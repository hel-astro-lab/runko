
#include "core/pic2/particle.h"

#include "tyvi/mdgrid.h"
#include "tyvi/mdspan.h"

namespace pic2 {

ParticleContainer::ParticleContainer(const ParticleContainerArgs args) :
  pos_(args.N),
  vel_(args.N),
  weights_(args.N),
  charge_ { args.charge },
  mass_ { args.mass }
{
}

std::size_t
  ParticleContainer::size() const
{
  return pos_.extents().extent(0);
}

ParticleContainerArgs
  ParticleContainer::args() const
{
  return { .N = size(), .charge = charge_, .mass = mass_ };
}

std::array<std::vector<ParticleContainer::value_type>, 3>
  ParticleContainer::get_positions()
{
  auto w = tyvi::mdgrid_work {}.sync_to_staging(pos_);

  auto x = std::vector<value_type>(this->size());
  auto y = std::vector<value_type>(this->size());
  auto z = std::vector<value_type>(this->size());

  w.wait();

  const auto mds = pos_.staging_mds();

  for(const auto [i]: tyvi::sstd::index_space(mds)) {
    x[i] = mds[i][0];
    y[i] = mds[i][1];
    z[i] = mds[i][2];
  }

  return { std::move(x), std::move(y), std::move(z) };
}

std::array<std::vector<ParticleContainer::value_type>, 3>
  ParticleContainer::get_velocities()
{

  auto w = tyvi::mdgrid_work {}.sync_to_staging(vel_);

  auto vx = std::vector<value_type>(this->size());
  auto vy = std::vector<value_type>(this->size());
  auto vz = std::vector<value_type>(this->size());

  w.wait();

  const auto mds = vel_.staging_mds();

  for(const auto [i]: tyvi::sstd::index_space(mds)) {
    vx[i] = mds[i][0];
    vy[i] = mds[i][1];
    vz[i] = mds[i][2];
  }

  return { std::move(vx), std::move(vy), std::move(vz) };
}

std::vector<ParticleContainer::value_type>
  ParticleContainer::get_weights()
{
  auto w = tyvi::mdgrid_work {}.sync_to_staging(weights_);

  auto weights = std::vector<value_type>(this->size());

  w.wait();

  const auto mds = weights_.staging_mds();

  for(const auto [i]: tyvi::sstd::index_space(mds)) { weights[i] = mds[i][]; }

  return weights;
}

void
  ParticleContainer::add_particles(ParticleContainer& other)
{
  if(&other == this) {
    throw std::runtime_error {
      "Trying to add particles to ParticleContainer from itself."
    };
  }

  const auto Nprev  = this->size();
  const auto Nother = other.size();
  const auto Ntotal = Nprev + Nother;

  auto new_pos     = VecGrid(Ntotal);
  auto new_vel     = VecGrid(Ntotal);
  auto new_weights = ScalarGrid(Ntotal);

  const auto new_pos_mds     = new_pos.mds();
  const auto new_vel_mds     = new_vel.mds();
  const auto new_weights_mds = new_weights.mds();

  const auto prev_pos_mds     = this->pos_.mds();
  const auto prev_vel_mds     = this->vel_.mds();
  const auto prev_weights_mds = this->weights_.mds();

  auto wA1 = tyvi::mdgrid_work {}.for_each_index(
    prev_pos_mds,
    [=](const auto idx, const auto tidx) {
      new_pos_mds[idx][tidx] = prev_pos_mds[idx][tidx];
      new_vel_mds[idx][tidx] = prev_vel_mds[idx][tidx];
    });
  auto wA2 = tyvi::mdgrid_work {}.for_each_index(prev_weights_mds, [=](const auto idx) {
    new_weights_mds[idx][] = prev_weights_mds[idx][];
  });

  const auto other_pos_mds     = other.pos_.mds();
  const auto other_vel_mds     = other.vel_.mds();
  const auto other_weights_mds = other.weights_.mds();

  const auto where_other_go = std::tuple { Nprev, Ntotal };

  const auto other_in_new_pos_mds     = std::submdspan(new_pos_mds, where_other_go);
  const auto other_in_new_vel_mds     = std::submdspan(new_vel_mds, where_other_go);
  const auto other_in_new_weights_mds = std::submdspan(new_weights_mds, where_other_go);

  auto wB1 = tyvi::mdgrid_work {}.for_each_index(
    other_in_new_pos_mds,
    [=](const auto idx, const auto tidx) {
      other_in_new_pos_mds[idx][tidx] = other_pos_mds[idx][tidx];
      other_in_new_vel_mds[idx][tidx] = other_vel_mds[idx][tidx];
    });
  auto wB2 =
    tyvi::mdgrid_work {}.for_each_index(other_in_new_weights_mds, [=](const auto idx) {
      other_in_new_weights_mds[idx][] = other_weights_mds[idx][];
    });

  tyvi::when_all(wA1, wA2, wB1, wB2).wait();

  pos_     = std::move(new_pos);
  vel_     = std::move(new_vel);
  weights_ = std::move(new_weights);
}
}  // namespace pic2

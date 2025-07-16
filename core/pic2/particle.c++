
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

}  // namespace pic2

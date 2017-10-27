
#include "sheets.h"



using namespace sheets;



/// Resize the sheet into correct size
void Sheet::resize(size_t Ni_, size_t Nj_) {
  iGrid.resize(Ni_);
  jGrid.resize(Nj_);
  values.resize(Ni_*Nj_);

  Ni = Ni_;
  Nj = Nj_;
}

/// internal function to get general id from sheet indices
size_t Sheet::getIndex(size_t i, size_t j) {
  return Ni*j + i;
}

/// Load scalar to the sheet
void Sheet::loadValue(size_t i, size_t j, Realf val) {
  size_t indx = getIndex(i, j);
  values[indx] = val;
}

void Sheet::loadZeroBlock(size_t i, size_t j) {
  size_t indx = getIndex(i, j);

  // FIXME add block instead of scalar
  values[indx] = 0.0;
}

void Sheet::loadBlock(size_t i, size_t j, vblock_t block) {
  size_t indx = getIndex(i, j);

  // FIXME add block instead of scalar
  values[indx] = block[0];
}

// return block at location (i,j)
vblock_t Sheet::getBlock(size_t i, size_t j) {
  vblock_t ret;
  size_t indx = getIndex(i, j);
  ret[0] = values[indx]; //FIXME return block instead element

  return ret;
}

bool Sheet::isNonZero(size_t i, size_t j) {
  size_t indx = getIndex(i, j);
  if ( values[indx] == 0.0 ) { return false; };
  return true;
};

// Sheet arithmetics
Sheet& Sheet::operator+=(const Sheet& rhs) {
  for(size_t q=0; q<(this->Ni*this->Nj); q++) this->values[q] += rhs.values[q];
  return *this;
}

Sheet& Sheet::operator-=(const Sheet& rhs) {
  for(size_t q=0; q<(this->Ni*this->Nj); q++) this->values[q] -= rhs.values[q];
  return *this;
}

Sheet& Sheet::operator*=(const Realf rhs) {
  for(size_t q=0; q<(this->Ni*this->Nj); q++) this->values[q] *= rhs;
  return *this;
}



/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Implementation for a Field.  Field is not intended so much to hide implementation
of data as to restrict write access to data.  It freely passes out pointers to
its private data, but only passes out read-only const pointers unless you have
the secret password (AKA the name of the process kernel that owns the data).

Field also stores some basic metadata for Vis, checkpointing, etc.
------------------------------------------------------------------------- */

#include "errors.hh"
#include "cell_geometry.hh"

#include "Field.hh"

namespace Amanzi {

Field::Field(std::string fieldname, FieldLocation location,
             Teuchos::RCP<Amanzi::AmanziMesh::Mesh> &mesh_maps,
             std::string owner, int num_dofs):
    fieldname_(fieldname), location_(location),
    owner_(owner), mesh_maps_(mesh_maps), num_dofs_(num_dofs),
    io_restart_(true), io_vis_(false) {

  if (location_ == FIELD_LOCATION_FACE) {
    data_ = Teuchos::rcp(new Epetra_MultiVector(mesh_maps->face_map(false),num_dofs_));
  } else if (location == FIELD_LOCATION_CELL) {
    data_ = Teuchos::rcp(new Epetra_MultiVector(mesh_maps->cell_map(false),num_dofs_));
  } else if (location == FIELD_LOCATION_MESH) {
    Errors::Message message("Constant Fields not yet implemented");
    Exceptions::amanzi_throw(message);
  }

  if (num_dofs == 1 && location == FIELD_LOCATION_CELL) {
    subfieldnames_.push_back(fieldname);
  }
};

// copy constructor:
// Create a new Field with different data but the same values.
Field::Field(const Field& field) {
  data_ = Teuchos::rcp(new Epetra_MultiVector(*field.data_));
  fieldname_ = field.fieldname_;
  owner_ = field.owner_;
  subfieldnames_ = field.subfieldnames_;
  location_ = field.location_;
  num_dofs_ = field.num_dofs_;
  io_restart_ = field.io_restart_;
  io_vis_ = field.io_vis_;
  initialized_ = field.initialized_;
  mesh_maps_ = field.mesh_maps_;
};

// assignment is crucial... copy values not data/pointers
Field& Field::operator=(const Field& field) {
  if (this != &field) {
    *data_ = *(field.data_);
    fieldname_ = field.fieldname_;
    owner_ = field.owner_;
    subfieldnames_ = field.subfieldnames_;
    location_ = field.location_;
    num_dofs_ = field.num_dofs_;
    io_restart_ = field.io_restart_;
    io_vis_ = field.io_vis_;
    initialized_ = field.initialized_;
    mesh_maps_ = field.mesh_maps_;
  }
  return *this;
};

// I miss decorators, do they exist in C++?
// check that the requesting pk owns the data
void Field::assert_owner_or_die(std::string pk_name) const {
  if (pk_name != owner_) {
    std::stringstream messagestream;
    messagestream << "PK " << pk_name << " is attempting to write to " << fieldname_
                  << " which is owned by " << owner_;
    Errors::Message message(messagestream.str());
    Exceptions::amanzi_throw(message);
  }
};

// check that the location matches this field's location
void Field::assert_location_or_die(FieldLocation location) const {
  if (location != location_) {
    std::stringstream messagestream;
    messagestream << "Invalid write to location " << location
                  << " when data is defined on " << location_;
    Errors::Message message(messagestream.str());
    Exceptions::amanzi_throw(message);
  }
};

// check that the block id is valid within the mesh
void Field::assert_valid_block_id_or_die(int mesh_block_id) const {
  // check valid block id
  if (!mesh_maps_->valid_set_id(mesh_block_id,Amanzi::AmanziMesh::CELL)) {
    std::stringstream messagestream;
    messagestream << "Invalid mesh block id: " << mesh_block_id;
    Errors::Message message(messagestream.str());
    Exceptions::amanzi_throw(message);
  }
};

// write-access to the data
Teuchos::RCP<Epetra_MultiVector> Field::get_data(std::string pk_name) {
  if (pk_name == owner_) {
    return data_;
  } else {
    std::stringstream messagestream;
    messagestream << "PK: " << pk_name << " is requesting write access to " << fieldname_
                  << " which is owned by " << owner_;
    Errors::Message message(messagestream.str());
    Exceptions::amanzi_throw(message);
  }
};

// write data
void Field::set_data(std::string pk_name, const Epetra_MultiVector& data) {
  assert_owner_or_die(pk_name);
  *data_ = data;
};

// write data
void Field::set_data(std::string pk_name, const Epetra_Vector& data) {
  assert_owner_or_die(pk_name);
  *((*data_)(0)) = data;
};

// write data to uniform, constant value (per vector)
void Field::set_data(std::string pk_name, const double* u) {
  assert_owner_or_die(pk_name);
  assert_location_or_die(FIELD_LOCATION_CELL);

  for (int i = 0; i != data_->NumVectors(); ++i) {
    ((*data_)(i))->PutScalar(u[i]);
  }
};

// write data to uniform, constant value
void Field::set_data(std::string pk_name, double u) {
  assert_owner_or_die(pk_name);
  assert_location_or_die(FIELD_LOCATION_CELL);

  for (int i = 0; i != data_->NumVectors(); ++i) {
    (*data_)(i)->PutScalar(u);
  }
};

// write data of on single block to a constant value (per vector)
void Field::set_data(std::string pk_name, const double* u, int mesh_block_id) {
  assert_owner_or_die(pk_name);
  assert_location_or_die(FIELD_LOCATION_CELL);
  assert_valid_block_id_or_die(mesh_block_id);

  unsigned int mesh_block_size = mesh_maps_->get_set_size(mesh_block_id,
                                                         Amanzi::AmanziMesh::CELL,
                                                         Amanzi::AmanziMesh::OWNED);
  std::vector<unsigned int> cell_ids(mesh_block_size);
  mesh_maps_->get_set(mesh_block_id, Amanzi::AmanziMesh::CELL, Amanzi::AmanziMesh::OWNED,
                     cell_ids.begin(),cell_ids.end());

  // loop over vectors and cells, assigning value
  for (int i = 0; i != data_->NumVectors(); ++i) {
    for( std::vector<unsigned int>::iterator c = cell_ids.begin();
         c != cell_ids.end();  c++) {
      (*(*data_)(i))[*c] = u[i];
    }
  }
};

// write data of on single block to a constant value
void Field::set_data(std::string pk_name, double u, int mesh_block_id) {
  assert_owner_or_die(pk_name);
  assert_location_or_die(FIELD_LOCATION_CELL);
  assert_valid_block_id_or_die(mesh_block_id);

  unsigned int mesh_block_size = mesh_maps_->get_set_size(mesh_block_id,
                                                         Amanzi::AmanziMesh::CELL,
                                                         Amanzi::AmanziMesh::OWNED);
  std::vector<unsigned int> cell_ids(mesh_block_size);
  mesh_maps_->get_set(mesh_block_id, Amanzi::AmanziMesh::CELL, Amanzi::AmanziMesh::OWNED,
                     cell_ids.begin(),cell_ids.end());

  // loop over vectors and cells, assigning value
  for (int i = 0; i != data_->NumVectors(); ++i) {
    for( std::vector<unsigned int>::iterator c = cell_ids.begin();
         c != cell_ids.end();  c++) {
      (*(*data_)(i))[*c] = u;
    }
  }
};

// write a constant 3-vector to faces (dotting the constant vector with the face normal)
void Field::set_vector_data(std::string pk_name, const double* u, int mesh_block_id) {
  assert_owner_or_die(pk_name);
  assert_location_or_die(FIELD_LOCATION_FACE);
  assert_valid_block_id_or_die(mesh_block_id);

  double x[4][3], normal[3];

  // get the cell ids
  unsigned int mesh_block_size = mesh_maps_->get_set_size(mesh_block_id,
                                                         Amanzi::AmanziMesh::CELL,
                                                         Amanzi::AmanziMesh::OWNED);
  std::vector<unsigned int> cell_ids(mesh_block_size);
  mesh_maps_->get_set(mesh_block_id, Amanzi::AmanziMesh::CELL, Amanzi::AmanziMesh::OWNED,
                     cell_ids.begin(),cell_ids.end());

  // loop over faces(cells), assigning value dotted with normal
  for( std::vector<unsigned int>::iterator c = cell_ids.begin();
       c != cell_ids.end();  c++) {

    std::vector<unsigned int> cface(6);
    mesh_maps_->cell_to_faces(*c, cface.begin(), cface.end());
    for (std::vector<unsigned int>::iterator f = cface.begin(); f != cface.end(); f++) {
      if (mesh_maps_->face_map(false).MyLID(*f) ) {
        mesh_maps_->face_to_coordinates( *f, (double*) x, (double*) x+12 );
        cell_geometry::quad_face_normal(x[0], x[1], x[2], x[3], normal);

        (*(*data_)(0))[*f] = u[0] * normal[0] + u[1] * normal[1] + u[2] * normal[2];
      }
    }
  }
};

} // namespace Amanzi

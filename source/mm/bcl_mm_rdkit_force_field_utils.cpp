// (c) Copyright BCL @ Vanderbilt University 2014
// (c) BCL Homepage: http://www.meilerlab.org/bclcommons
// (c) BCL Code Repository: https://github.com/BCLCommons/bcl
// (c)
// (c) The BioChemical Library (BCL) was originally developed by contributing members of the Meiler Lab @ Vanderbilt University.
// (c)
// (c) The BCL is now made available as an open-source software package distributed under the permissive MIT license,
// (c) developed and maintained by the Meiler Lab at Vanderbilt University and contributing members of the BCL Commons.
// (c)
// (c) External code contributions to the BCL are welcome. Please visit the BCL Commons GitHub page for information on how you can contribute.
// (c)
// (c) This file is part of the BCL software suite and is made available under the MIT license.
// (c)

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "mm/bcl_mm_rdkit_force_field_utils.h"

// includes from bcl - sorted alphabetically
#include "storage/bcl_storage_map.h"
#include "util/bcl_util_message.h"
#include "util/bcl_util_undefined.h"

// external includes - sorted alphabetically
#include "ForceField/MMFF/PositionConstraint.h"
#include "ForceField/UFF/PositionConstraint.h"
#include "GraphMol/ForceFieldHelpers/MMFF/MMFF.h"
#include "GraphMol/ForceFieldHelpers/MMFF/AtomTyper.h"
#include "GraphMol/ForceFieldHelpers/MMFF/Builder.h"
#include "GraphMol/ForceFieldHelpers/UFF/AtomTyper.h"
#include "GraphMol/ForceFieldHelpers/UFF/Builder.h"

namespace bcl
{
  namespace mm
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief virtual copy constructor
    RdkitForceFieldUtils *RdkitForceFieldUtils::Clone() const
    {
      return new RdkitForceFieldUtils( *this);
    }

    //! @brief returns class name
    //! the class name as const ref std::string
    const std::string &RdkitForceFieldUtils::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

  ///////////////
  // operators //
  ///////////////

    //! @brief return the enum corresponding to the force field name
    //! @param NAME the string corresponding to the force field enum
    //! @return the force field enum
    RdkitForceFieldsEnum RdkitForceFieldUtils::GetRdkitForceFieldsEnum( const std::string &NAME)
    {
      static const std::map<std::string, RdkitForceFields> ff_map
      {
        { "UFF", e_UFF },
        { "MMFF94", e_MMFF94 },
        { "MMFF94s", e_MMFF94s }
      };
      auto itr( ff_map.find( NAME));
      return ( itr != ff_map.end() ) ? itr->second : s_NumberForceFields;
    }

    //! @brief return an MMFF force field
    //! @param RDMOL the rdkit mol to be parameterized
    //! @param MMFF_VARIANT specify MMFF94 or MMFF94s
    //! @param NONBONDED_THRESHOLD the threshold to be used in adding non-bonded terms to the force field.
    //! Any non-bonded contact whose current distance is greater than nonBondedThresh * the minimum value for that contact
    //! will not be included.
    //! @param IGNORE_INTER_FRAGMENT_INTERACTIONS if true, nonbonded terms will not be added between fragments
    //! @param INITIALIZE if true, initialize the force field for use
    //! @return a force field
    ::ForceFields::ForceField* RdkitForceFieldUtils::ConstructForceField
    (
      ::RDKit::RWMol &RDMOL,
      const std::string &MMFF_VARIANT,
      const double NONBONDED_THRESHOLD,
      const bool IGNORE_INTER_FRAGMENT_INTERACTIONS,
      const bool INITIALIZE
    )
    {
      if( GetRdkitForceFieldsEnum( MMFF_VARIANT) == e_UFF)
      {
        ::ForceFields::ForceField *ff = ::RDKit::UFF::constructForceField( RDMOL, NONBONDED_THRESHOLD, -1, IGNORE_INTER_FRAGMENT_INTERACTIONS);
        if( INITIALIZE){ ff->initialize(); }
        return ff;
      }
      else if( GetRdkitForceFieldsEnum( MMFF_VARIANT) == e_MMFF94 || GetRdkitForceFieldsEnum( MMFF_VARIANT) == e_MMFF94s)
      {
        ::RDKit::MMFF::MMFFMolProperties mmff_mol_properties( RDMOL, MMFF_VARIANT);
        BCL_Assert( mmff_mol_properties.isValid(), "Invalid MMFF molecule properties!");
        ::ForceFields::ForceField *ff = ::RDKit::MMFF::constructForceField( RDMOL, NONBONDED_THRESHOLD, -1, IGNORE_INTER_FRAGMENT_INTERACTIONS);
        if( INITIALIZE){ ff->initialize(); }
        return ff;
      }

      // add more here as needed / implemented

      // defaulting
      BCL_MessageStd("Warning: no valid RDKitForceField specified; defaulting to UFF");
      ::ForceFields::ForceField *ff = ::RDKit::UFF::constructForceField( RDMOL, NONBONDED_THRESHOLD, -1, IGNORE_INTER_FRAGMENT_INTERACTIONS);
      if( INITIALIZE){ ff->initialize(); }
      return ff;
    }

    //! @brief add positional constraints to force field for geometry optimization
    //! @param FORCE_FIELD the force field that is modified with the new restraint term
    //! @param ATOM_INDICES indices that are restrained during minimization
    //! @param MAX_UNRESTRAINED DISPLACEMENT coordinate displacement above which restraint force is applied
    //! @param RESTRAINT_FORCE restraint force
    void RdkitForceFieldUtils::AddPositionalRestraintsUFF
    (
      ::ForceFields::ForceField *FORCE_FIELD, // raw pointer unconventional for BCL outside of Clone(), but this is what RDKit requires
      const storage::Vector< size_t> &ATOM_INDICES,
      const storage::Vector< double> &MAX_UNRESTRAINED_DISPLACEMENT,
      const storage::Vector< double> &RESTRAINT_FORCE
    )
    {
      // sanity check on vector sizes
      if( ATOM_INDICES.GetSize() == MAX_UNRESTRAINED_DISPLACEMENT.GetSize() && ATOM_INDICES.GetSize() == RESTRAINT_FORCE.GetSize())
      {
        // loop over atom indices and add restraint forces to our force field
        for( size_t i( 0), sz( ATOM_INDICES.GetSize()); i < sz; ++i)
        {
          ::ForceFields::UFF::PositionConstraintContrib *coord_cst;
          coord_cst = new ::ForceFields::UFF::PositionConstraintContrib
              (
                FORCE_FIELD, ATOM_INDICES( i),
                MAX_UNRESTRAINED_DISPLACEMENT( i),
                RESTRAINT_FORCE( i)
              );
          FORCE_FIELD->contribs().push_back( ForceFields::ContribPtr( coord_cst));
        }
      }
      else if( ATOM_INDICES.GetSize() == MAX_UNRESTRAINED_DISPLACEMENT.GetSize() && RESTRAINT_FORCE.GetSize() == size_t( 1))
      {
        for( size_t i( 0), sz( ATOM_INDICES.GetSize()); i < sz; ++i)
        {
          ::ForceFields::UFF::PositionConstraintContrib *coord_cst;
          coord_cst = new ::ForceFields::UFF::PositionConstraintContrib
              (
                FORCE_FIELD, ATOM_INDICES( i),
                MAX_UNRESTRAINED_DISPLACEMENT( i),
                RESTRAINT_FORCE( 0)
              );
          FORCE_FIELD->contribs().push_back( ForceFields::ContribPtr( coord_cst));
        }
      }
      // do not kill, but inform user that positional restraints are not added
      else
      {
        BCL_MessageStd
        (
          "[WARNING] RdkitForceFieldUtils::AddPositionalRestraintsUFF "
          "The number of atoms does not match the number of max displacements and/or the number of provided restraint forces; "
          "alternatively, if the number of restraint forces to be added is one, then the number of atoms simply does not match the number of "
          "max displacements. NO POSITIONAL RESTRAINT ADDED!"
          " Number of atoms to be restrained: " + util::Format()( ATOM_INDICES.GetSize()) +
          " Number of max unrestrained distances: " + util::Format()( MAX_UNRESTRAINED_DISPLACEMENT.GetSize()) +
          " Number of restraint forces specified: " + util::Format()( RESTRAINT_FORCE.GetSize())
        )
      }
    }

    //! @brief add positional constraints to force field for geometry optimization
    //! @param FORCE_FIELD the force field that is modified with the new restraint term
    //! @param ATOM_INDICES indices that are restrained during minimization
    //! @param MAX_UNRESTRAINED DISPLACEMENT coordinate displacement above which restraint force is applied
    //! @param RESTRAINT_FORCE restraint force
    void RdkitForceFieldUtils::AddPositionalRestraintsMMFF
    (
      ::ForceFields::ForceField *FORCE_FIELD, // raw pointer unconventional for BCL outside of Clone(), but this is what RDKit requires
      const storage::Vector< size_t> &ATOM_INDICES,
      const storage::Vector< double> &MAX_UNRESTRAINED_DISPLACEMENT,
      const storage::Vector< double> &RESTRAINT_FORCE
    )
    {
      // sanity check on vector sizes
      if( ATOM_INDICES.GetSize() == MAX_UNRESTRAINED_DISPLACEMENT.GetSize() && ATOM_INDICES.GetSize() == RESTRAINT_FORCE.GetSize())
      {
        // loop over atom indices and add restraint forces to our force field
        for( size_t i( 0), sz( ATOM_INDICES.GetSize()); i < sz; ++i)
        {
          ::ForceFields::MMFF::PositionConstraintContrib *coord_cst;
          coord_cst = new ::ForceFields::MMFF::PositionConstraintContrib
              (
                FORCE_FIELD, ATOM_INDICES( i),
                MAX_UNRESTRAINED_DISPLACEMENT( i),
                RESTRAINT_FORCE( i)
              );
          FORCE_FIELD->contribs().push_back( ForceFields::ContribPtr( coord_cst));
        }
      }
      else if( ATOM_INDICES.GetSize() == MAX_UNRESTRAINED_DISPLACEMENT.GetSize() && RESTRAINT_FORCE.GetSize() == size_t( 1))
      {
        for( size_t i( 0), sz( ATOM_INDICES.GetSize()); i < sz; ++i)
        {
          ::ForceFields::MMFF::PositionConstraintContrib *coord_cst;
          coord_cst = new ::ForceFields::MMFF::PositionConstraintContrib
              (
                FORCE_FIELD, ATOM_INDICES( i),
                MAX_UNRESTRAINED_DISPLACEMENT( i),
                RESTRAINT_FORCE( 0)
              );
          FORCE_FIELD->contribs().push_back( ForceFields::ContribPtr( coord_cst));
        }
      }
      // do not kill, but inform user that positional restraints are not added
      else
      {
        BCL_MessageStd
        (
          "[WARNING] RdkitForceFieldUtils::AddPositionalRestraintsMMFF "
          "The number of atoms does not match the number of max displacements and/or the number of provided restraint forces; "
          "alternatively, if the number of restraint forces to be added is one, then the number of atoms simply does not match the number of "
          "max displacements. NO POSITIONAL RESTRAINT ADDED!"
          " Number of atoms to be restrained: " + util::Format()( ATOM_INDICES.GetSize()) +
          " Number of max unrestrained distances: " + util::Format()( MAX_UNRESTRAINED_DISPLACEMENT.GetSize()) +
          " Number of restraint forces specified: " + util::Format()( RESTRAINT_FORCE.GetSize())
        )
      }
    }


  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &RdkitForceFieldUtils::Read( std::istream &ISTREAM)
    {
      // read member
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return output stream which was written to
    std::ostream &RdkitForceFieldUtils::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

  } // namespace mm
} // namespace bcl


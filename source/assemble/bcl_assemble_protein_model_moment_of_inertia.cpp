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
#include "assemble/bcl_assemble_protein_model_moment_of_inertia.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_aa_neighbor_list_container.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "biol/bcl_biol_atom.h"
#include "coord/bcl_coord_moment_of_inertia.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> ProteinModelMomentOfInertia::s_Instance
    (
      GetObjectInstances().AddInstance
      (
        new ProteinModelMomentOfInertia
        (
          biol::AATypeData::e_TransferFreeEnergyEisenberg,
          util::ShPtr< AAExposureInterface>(),
          false
        )
      )
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    //! @param PROPERTY the property to use as weight
    //! @param EXPOSURE_WEIGHT_FUNCTION multiply each property with aa exposure
    //! @param CONSIDER_ONLY_NEGATIVE_PROPERTY consider only negative energies in calculations
    ProteinModelMomentOfInertia::ProteinModelMomentOfInertia
    (
      const biol::AATypeData::PropertyType &PROPERTY,
      const util::ShPtr< AAExposureInterface> &EXPOSURE_WEIGHT_FUNCTION,
      const bool CONSIDER_ONLY_NEGATIVE_PROPERTY
    ) :
      m_PropertyType( PROPERTY),
      m_ExposureWeight( EXPOSURE_WEIGHT_FUNCTION),
      m_ConsiderOnlyNegativeEnergies( CONSIDER_ONLY_NEGATIVE_PROPERTY)
    {
    }

    //! @brief Clone function
    //! @return pointer to new MomentOfInertia
    ProteinModelMomentOfInertia *ProteinModelMomentOfInertia::Clone() const
    {
      return new ProteinModelMomentOfInertia( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ProteinModelMomentOfInertia::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief return a 4*n matrix with 3 coordinates and exposure * weight in each row (4 columns) and number amino acid rows
    //! @param PROTEIN_MODEL protein model
    //! @return 4 * n matrix
    linal::Matrix< double>
    ProteinModelMomentOfInertia::ProteinModelAAExposureToCoordinateWeightMatrix( const ProteinModel &PROTEIN_MODEL) const
    {
      // collect AANeighborLists for all amino acids in the given protein model
      const AANeighborListContainer all_aa_neighbor_list
      (
        PROTEIN_MODEL.GetAminoAcids(),
        m_ExposureWeight->GetDistanceCutoff(),
        0,
        true
      );

      // create matrix
      linal::Matrix< double> coordinate_weight_matrix( all_aa_neighbor_list.GetSize(), 4, 0.0);

      // pointer to first row
      double *row( coordinate_weight_matrix.Begin());

      // iterate over all lists in the all_aa_neighbor_list
      for
      (
        AANeighborListContainer::const_iterator
          aa_itr( all_aa_neighbor_list.Begin()), aa_itr_end( all_aa_neighbor_list.End());
        aa_itr != aa_itr_end;
        ++aa_itr, row += 4
      )
      {
        const linal::Vector3D &current_coord( aa_itr->second.GetCenterAminoAcid()->GetFirstSidechainAtom().GetCoordinates());
        if( !current_coord.IsDefined())
        {
          continue;
        }

        // coordinate
        std::copy( current_coord.Begin(), current_coord.End(), row);

        // weight
        double exposure( m_ExposureWeight->operator ()( aa_itr->second));
        if( !m_ExposureWeight->IsDirect())
        {
          exposure = 1.0 / exposure;
        }
        const double current_weight( aa_itr->second.GetCenterAminoAcid()->GetType()->GetAAProperty( m_PropertyType) * exposure);
        if( !m_ConsiderOnlyNegativeEnergies || current_weight < 0.0)
        {
          row[ 3] = current_weight;
        }
      }

      // end
      return coordinate_weight_matrix;
    }

    //! @brief return a 4*n matrix with 3 coordinates and weight in each row (4 columns) and number amino acid rows
    //! @param PROTEIN_MODEL protein model
    //! @return 4 * n matrix
    linal::Matrix< double>
    ProteinModelMomentOfInertia::ProteinModelToCoordinateWeightMatrix( const ProteinModel &PROTEIN_MODEL) const
    {
      // get all amino acids of the protein model
      const util::SiPtrVector< const biol::AABase> amino_acids( PROTEIN_MODEL.GetAminoAcids());

      // create matrix
      linal::Matrix< double> coordinate_weight_matrix( amino_acids.GetSize(), 4, 0.0);

      // pointer to first row
      double *row( coordinate_weight_matrix.Begin());

      // iterate over all amino acids
      for
      (
        util::SiPtrVector< const biol::AABase>::const_iterator itr( amino_acids.Begin()), itr_end( amino_acids.End());
        itr != itr_end;
        ++itr, row += 4
      )
      {
        const linal::Vector3D &current_coord( ( *itr)->GetFirstSidechainAtom().GetCoordinates());
        if( !current_coord.IsDefined())
        {
          continue;
        }

        // coordinate
        std::copy( current_coord.Begin(), current_coord.End(), row);

        // weight
        const double current_weight( ( *itr)->GetType()->GetAAProperty( m_PropertyType));
        if( !m_ConsiderOnlyNegativeEnergies || current_weight < 0.0)
        {
          row[ 3] = current_weight;
        }
      }

      // end
      return coordinate_weight_matrix;
    }

    //! @brief calculate transformation that translates into the center of weights and that sorts principal axes of inertia according to principal moments of inertia x - smallest, z - largest
    //! @param PROTEIN_MODEL protein model
    //! @return transformation matrix and principal moments of inertias
    storage::Pair< math::TransformationMatrix3D, linal::Vector3D>
    ProteinModelMomentOfInertia::TransformationAndMoments( const ProteinModel &PROTEIN_MODEL) const
    {
      if( m_ExposureWeight.IsDefined())
      {
        return coord::MomentOfInertia().TransformationAndMoments( ProteinModelAAExposureToCoordinateWeightMatrix( PROTEIN_MODEL));
      }
      else
      {
        return coord::MomentOfInertia().TransformationAndMoments( ProteinModelToCoordinateWeightMatrix( PROTEIN_MODEL));
      }
    }

  ///////////////
  // operators //
  ///////////////

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ProteinModelMomentOfInertia::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_PropertyType, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &ProteinModelMomentOfInertia::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_PropertyType, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  } // namespace assemble
} // namespace bcl

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
#include "restraint/bcl_restraint_rdc.h"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_matrix3x3.h"
#include "linal/bcl_linal_matrix_operations.h"
#include "linal/bcl_linal_vector_3d_operations.h"
#include "restraint/bcl_restraint_rdc_assignment.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> RDC::s_Instance
    (
      GetObjectInstances().AddInstance( new RDC())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    RDC::RDC() :
      m_Data()
    {
    }

    //! @brief Clone function
    //! @return pointer to new RDC
    RDC *RDC::Clone() const
    {
      return new RDC( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &RDC::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief normalizes the RDC value to an NH value
    void RDC::NormalizetoNH()
    {
      // calculate the NH normalization factor
      static const double s_nh_value
      (
        CalculateNormalization
        (
          chemistry::GetElementTypes().e_Nitrogen,
          chemistry::GetElementTypes().e_Hydrogen,
          biol::GetAtomTypes().N->GetBondLength( biol::GetAtomTypes().H)
        )
      );

      // iterate through the rdc data
      for
      (
        storage::Vector< storage::Triplet< DataPairwise, double, double> >::iterator
          itr( m_Data.Begin()), itr_end( m_Data.End()); itr != itr_end; ++itr
      )
      {
        // update the rdc value
        itr->Third() /=
          CalculateNormalization
          (
            itr->First().First()->GetAtomType()->GetElementType(),  // first element type
            itr->First().Second()->GetAtomType()->GetElementType(), // second element type
            itr->Second()                                        // internuclear distance
          ) /
          s_nh_value;
      }
    }

    //! @brief adjust the RDC signs relative to NH (i.e. value is correct, but sign is flipped)
    void RDC::AdjustSigns()
    {
      // iterate through the rdc data
      for
      (
        storage::Vector< storage::Triplet< DataPairwise, double, double> >::iterator
          itr( m_Data.Begin()), itr_end( m_Data.End()); itr != itr_end; ++itr
      )
      {
        // if the two gyromagnetic ratios are positive when multiplied (NH will be negative)
        if
        (
          itr->First().First()->GetAtomType()->GetElementType()->GetProperty
          (
            chemistry::ElementTypeData::e_GyromagneticRatio
          ) *
          itr->First().Second()->GetAtomType()->GetElementType()->GetProperty
          (
            chemistry::ElementTypeData::e_GyromagneticRatio
          ) >
          0.0
        )
        {
          // switch the sign of the rdc value
          itr->Third() = -itr->Third();
        }
      }
    }

    //! @brief calculates an RDC value from two coordinates and a tensor
    //! @param COORDINATES_A first coordinates
    //! @param COORDINATES_B second coordinates
    //! @param TENSOR tensor to be applied
    //! @return calculated RDC value
    double RDC::CalculateValue
    (
      const linal::Vector3D &COORDINATES_A,
      const linal::Vector3D &COORDINATES_B,
      const linal::Matrix3x3< double> &TENSOR
    )
    {
      // calculate the vector
      linal::Vector3D vector_value( linal::UnitVector( COORDINATES_B, COORDINATES_A));

      // store the values as a matrix
      const linal::Matrix< double> coordinates( 3, 1, vector_value.Begin());

      // multiply the transposed coordinates by the tensor and then by themselves again in order
      // to obtain an RDC value
      return ( coordinates.Transposed() * TENSOR * coordinates)( 0, 0);
    }

    //! @brief add new data into the Vector
    //! @param LOCATOR_A first atom locator
    //! @param LOCATOR_B second atom locator
    //! @param INTERNUCLEAR_DISTANCE internuclear distance
    //! @param VALUE experimental RDC value
    void RDC::PushBack
    (
      const assemble::LocatorAtomCoordinatesInterface &LOCATOR_A,
      const assemble::LocatorAtomCoordinatesInterface &LOCATOR_B,
      const double INTERNUCLEAR_DISTANCE,
      const double VALUE
    )
    {
      m_Data.PushBack
      (
        storage::Triplet< DataPairwise, double, double>
        (
          DataPairwise
          (
            util::ShPtr< assemble::LocatorAtomCoordinatesInterface>( LOCATOR_A.Clone()),
            util::ShPtr< assemble::LocatorAtomCoordinatesInterface>( LOCATOR_B.Clone())
          ),
          INTERNUCLEAR_DISTANCE,
          VALUE
        )
      );
    }

    //! @brief shuffles the Vector of RDCs
    void RDC::Shuffle()
    {
      std::random_shuffle( m_Data.Begin(), m_Data.End());
    }

    //! @brief generates the assignment from the protein model
    //! @param PROTEIN_MODEL protein model to be used to generate the assignment
    //! @return assignment containing coordinates of located atoms and experimental distance
    RDCAssignment RDC::GenerateAssignment( const assemble::ProteinModel &PROTEIN_MODEL) const
    {
      // initialize RDC assignment
      RDCAssignment assignment;

      // iterate through the rdc data
      for
      (
        storage::Vector< storage::Triplet< DataPairwise, double, double> >::const_iterator
          itr( m_Data.Begin()), itr_end( m_Data.End()); itr != itr_end; ++itr
      )
      {
        // get the atom coordinates
        const linal::Vector3D first_atom_coords( itr->First().First()->Locate( PROTEIN_MODEL));
        const linal::Vector3D second_atom_coords( itr->First().Second()->Locate( PROTEIN_MODEL));
        const double value( itr->Third());

        // if all the values are defined and the experimental value is non-zero
        if( first_atom_coords.IsDefined() && second_atom_coords.IsDefined() && util::IsDefined( value) && value != 0.0)
        {
          // add the information to the assignment
          assignment.PushBack( first_atom_coords, second_atom_coords, value);
        }
      }

      // end
      return assignment;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &RDC::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Data, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &RDC::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Data, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief calculate the normalization factor for the given element types and bond length
    //! @param ELEMENT_TYPE_A first element type
    //! @param ELEMENT_TYPE_B first element type
    //! @param BOND_LENGTH bond length
    //! @return the normalization factor
    double RDC::CalculateNormalization
    (
      const chemistry::ElementType &ELEMENT_TYPE_A,
      const chemistry::ElementType &ELEMENT_TYPE_B,
      const double BOND_LENGTH
    )
    {
      // calculate the normalization factor
      return ELEMENT_TYPE_A->GetProperty( chemistry::ElementTypeData::e_GyromagneticRatio) *
        ELEMENT_TYPE_B->GetProperty( chemistry::ElementTypeData::e_GyromagneticRatio) /
        math::Pow( BOND_LENGTH, 3.0);
    }

  } // namespace restraint

} // namespace bcl

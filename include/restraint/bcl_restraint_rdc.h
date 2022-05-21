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

#ifndef BCL_RESTRAINT_RDC_H_
#define BCL_RESTRAINT_RDC_H_

// include the namespace header
#include "bcl_restraint.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_restraint_data_pairwise.h"
#include "assemble/bcl_assemble_locator_atom.h"
#include "command/bcl_command_flag_interface.h"
#include "storage/bcl_storage_list.h"
#include "storage/bcl_storage_triplet.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class RDC
    //! @brief NMR RDC restraint class
    //! @details Contains RDC restraint information for a protein model in the form of two LocatorAtoms, an
    //!          internuclear distance, and an experimental value.  A collections of RDCs are stored since each
    //!          restraint is only meaningful in the context of all the other restraints.  Functions are included to
    //!          normalize the data based on the element types of the atoms.
    //!
    //! @see @link example_restraint_rdc.cpp @endlink
    //! @author weinerbe
    //! @date Feb 16, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API RDC :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      //! list of RDC restraints, 2 atom locators, the internuclear distance, and the experimental value
      storage::Vector< storage::Triplet< DataPairwise, double, double> > m_Data;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      RDC();

      //! @brief Clone function
      //! @return pointer to new RDC
      RDC *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief get the rdc data
      //! @return the rdc data
      const storage::Vector< storage::Triplet< DataPairwise, double, double> > &GetData() const
      {
        return m_Data;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief normalizes the RDC value to an NH value
      void NormalizetoNH();

      //! @brief adjust the RDC sign relative to NH (i.e. value is correct, but sign is flipped)
      void AdjustSigns();

      //! @brief add new data into the list
      //! @param LOCATOR_A first atom locator
      //! @param LOCATOR_B second atom locator
      //! @param INTERNUCLEAR_DISTANCE internuclear distance
      //! @param VALUE experimental RDC value
      void PushBack
      (
        const assemble::LocatorAtomCoordinatesInterface &LOCATOR_A,
        const assemble::LocatorAtomCoordinatesInterface &LOCATOR_B,
        const double INTERNUCLEAR_DISTANCE,
        const double VALUE
      );

      //! @brief shuffles the Vector of RDCs
      void Shuffle();

      //! @brief generates the assignment from the protein model
      //! @param PROTEIN_MODEL protein model to be used to generate the assignment
      //! @return assignment containing coordinates of located atoms and experimental distance
      RDCAssignment GenerateAssignment( const assemble::ProteinModel &PROTEIN_MODEL) const;

      //! @brief calculates an RDC value from two coordinates and a tensor
      //! @param COORDINATES_A first coordinates
      //! @param COORDINATES_B second coordinates
      //! @param TENSOR tensor to be applied
      //! @return calculated RDC value
      static double CalculateValue
      (
        const linal::Vector3D &COORDINATES_A,
        const linal::Vector3D &COORDINATES_B,
        const linal::Matrix3x3< double> &TENSOR
      );

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

    private:

      //! @brief calculate the normalization factor for the given element types and bond length
      //! @param ELEMENT_TYPE_A first element type
      //! @param ELEMENT_TYPE_B first element type
      //! @param BOND_LENGTH bond length
      //! @return the normalization factor
      static double CalculateNormalization
      (
        const chemistry::ElementType &ELEMENT_TYPE_A,
        const chemistry::ElementType &ELEMENT_TYPE_B,
        const double BOND_LENGTH
      );

    }; // class RDC

  } // namespace restraint
} // namespace bcl

#endif // BCL_RESTRAINT_RDC_H_

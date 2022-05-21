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

#ifndef BCL_RESTRAINT_ACCESSIBILITY_AA_H_
#define BCL_RESTRAINT_ACCESSIBILITY_AA_H_

// include the namespace header
#include "bcl_restraint.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"
#include "biol/bcl_biol.fwd.hh"
#include "storage/bcl_storage.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_aa_exposure_interface.h"
#include "assemble/bcl_assemble_locator_aa.h"
#include "math/bcl_math_range.h"
#include "storage/bcl_storage_set.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AccessibilityAA
    //! @brief AccessibilityAA is a restraint class which uses accessibility as the restraining data.
    //!        It can be used to indicate the accessibility that an amino acid should have according to the data.
    //!
    //! @see @link example_restraint_accessibility_aa.cpp @endlink
    //!
    //! @author alexanns
    //! @date 11/02/2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    class BCL_API AccessibilityAA :
      public util::ObjectInterface
    {
    public:

    /////////////////
    // environment //
    /////////////////

      enum EnvironmentType
      {
        e_Oxygen,
        e_NiEDDA,
        s_NumberEnvironmentTypes
      };

      //! @brief EnvironmentType as string
      //! @param TYPE the type
      //! @return the string for TYPE
      static const std::string &GetEnvironmentName( const EnvironmentType &TYPE);

      //! @brief EnvironmentType enum I/O helper
      typedef util::WrapperEnum< EnvironmentType, &GetEnvironmentName, s_NumberEnvironmentTypes> EnvironmentEnum;

    private:

    //////////
    // data //
    //////////

      //! the accessibility measured
      storage::Map< EnvironmentEnum, double> m_Accessibility;

      //! the amino acid the accessibility corresponds to
      util::ShPtr< assemble::LocatorAA> m_AminoAcid;

      //! the method to use for calculating exposure of residues from structure
      util::ShPtr< assemble::AAExposureInterface> m_ExposureCalculator;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      AccessibilityAA();

      //! @brief construct from all the member variables
      //! @param ACCESSIBILITY double which is the accessibility
      //! @param AMINO_ACID the type of accessibility measurement made
      //! @param EXPOSURE_CALCULATOR the method to use for calculating exposure of residues from structure
      AccessibilityAA
      (
        const storage::Map< EnvironmentEnum, double> &ACCESSIBILITY,
        const util::ShPtr< assemble::LocatorAA> &AMINO_ACID,
        const util::ShPtr< assemble::AAExposureInterface> &EXPOSURE_CALCULATOR
      );

      //! @brief virtual copy constructor
      virtual AccessibilityAA *Clone() const;

      //! @brief virtual destructor
      virtual ~AccessibilityAA();

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      virtual const std::string &GetClassIdentifier() const;

      //! @brief gives a const reference to the residue locator
      //! @return a const reference to the residue locator
      const util::ShPtr< assemble::LocatorAA> &GetAA() const
      {
        return m_AminoAcid;
      }

      //! @brief gives const reference to the set of accessibilities
      //! @return const reference to the set of accessibilities
      const storage::Map< EnvironmentEnum, double> &GetAccessibilityAAs() const
      {
        return m_Accessibility;
      }

      //! @brief provides the desired accessibility based on a given environment type
      //! @return pair of bool and Accessibility where the bool indicates if that environment type exists and
      //!         the Accessibility is the experimentally measured value
      storage::Pair< bool, double>
      GetAccessibilityByEnvironment( const AccessibilityAA::EnvironmentType &ENVIRONMENT) const;

    ////////////////
    // operations //
    ////////////////

      //! @brief GenerateNeighborList generates an AANeighborList from a protein model for the desired residue
      //! @param PROTEIN_MODEL is used to create the assemble::AANeighborList
      //! @return returns assemble::AANeighborList for the located residue from protein model
      assemble::AANeighborList
      GenerateNeighborList
      (
        const assemble::ProteinModel &PROTEIN_MODEL
      ) const;

      //! @brief GenerateAssignment generates an assignment for the located residue from a protein model
      //! @param PROTEIN_MODEL is used to create the AccessibilityAAAssignment
      //! @return returns AccessibilityAAAssignment for the located residue from protein model
      AccessibilityAAAssignment GenerateAssignment
      (
        const assemble::ProteinModel &PROTEIN_MODEL
      ) const;

    //////////////////////
    // input and output //
    //////////////////////
    protected:

      //! @brief read distance from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      virtual std::istream &Read( std::istream &ISTREAM);

      //! @brief write distance to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT indentation
      //! @return output stream which was written to
      virtual std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; //< AccessibilityAA

  } // namespace restraint
} // namespace bcl

#endif //BCL_RESTRAINT_ACCESSIBILITY_AA_H_

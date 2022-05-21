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

#ifndef BCL_APP_GENERATE_ATOM_HYBRIDIZATION_DESCRIPTORS_H_
#define BCL_APP_GENERATE_ATOM_HYBRIDIZATION_DESCRIPTORS_H_

// include the interface for all apps
#include "app/bcl_app_apps.h"
#include "chemistry/bcl_chemistry_fragment_ensemble.h"

// include other forward headers - sorted alphabetically
#include "chemistry/bcl_chemistry.fwd.hh"
#include "math/bcl_math.fwd.hh"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_fragment_complete.h"
#include "command/bcl_command_command.h"
#include "command/bcl_command_flag_interface.h"
#include "linal/bcl_linal_vector.h"
#include "math/bcl_math_running_average.h"
#include "math/bcl_math_running_average_sd.h"
#include "storage/bcl_storage_vector.h"
#include "storage/bcl_storage_vector_nd.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace app
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class GenerateAtomHybridizationDescriptors
    //! @brief Application for generating atom descriptors to determine the hybridization of different atom types and
    //!        determining covalent and van der waals radii from a molecule structure library, such as the cambridge
    //!        structural database
    //!
    //! @author mendenjl, brownbp1
    //! @date Aug 1, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API GenerateAtomHybridizationDescriptors :
      public Interface
    {
    private:

    //////////
    // data //
    //////////

      //! Output filename base
      util::ShPtr< command::FlagInterface> m_OutputFilenameBase;

      //! output initializer for chemistry::BondLengths
      util::ShPtr< command::FlagInterface> m_WriteBondLengthsInitializer;

      //! output initializer for chemistry::BondLengths
      util::ShPtr< command::FlagInterface> m_WriteDescriptors;

      //! count element_type-element_type-bond_type statistics (normalized)
      util::ShPtr< command::FlagInterface> m_WriteElementBondStatistics;

      //! whether to skip aromatic rings for descriptors
      util::ShPtr< command::FlagInterface> m_SkipAromaticRingDescriptors;

      //! bond length statistics, keyed by ring size - 3 (outer vector), sdf bond type (vector) and atom types (inner map)
      mutable
        storage::Vector
        <
          storage::Vector< storage::Map< std::string, math::RunningAverageSD< double> > >
        > m_BondLengths;

      //! # of basic bond types
      static const size_t s_NumberBondTypes = 5;

      //! # of ring types; set this > 1 to split out differences in bond length by smallest-ring size membership
      //! Generally, this has no effect for C, for N and O the effect may be significant
      static const size_t s_NumberRingTypes = 1;

      //! All atom types using each basic bond type (key)
      //! Basic bond types:
      //! 1  - Single bond
      //! 2  - Double bond
      //! 3  - Triple bond
      //! 4  - Aromatic bond
      mutable storage::Vector< storage::Vector< storage::Set< std::string> > > m_AtomTypesForSimpleBondType;

      //! Datasets, indexed by string containing ElementSymbol # bonds, and # e- in bonds
      mutable storage::Map< std::string, storage::Vector< storage::VectorND< 2, linal::Vector< double> > > >
        m_AngleDescriptorsByType;

      mutable
        storage::Map< std::string, storage::Vector< double> > m_VdwBondLengths;
      mutable
        storage::Vector< storage::Map< std::string, storage::Vector< double> > > m_VdwBondLengthsByBondDist;
       mutable
        storage::Map< std::string, double> m_ExpectedLengthsByBondDist;
       mutable
        storage::Map< std::string, double> m_ExpectedVdw;
       mutable
        storage::Vector< storage::Map< std::string, double> > m_TimesNotTheClosest;

    ///////////////////////////////////
    // construction and destruction //
    ///////////////////////////////////

      //! default constructor
      GenerateAtomHybridizationDescriptors();

    public:

      // instantiate enumerator for GenerateAtomHybridizationDescriptors class
      static const ApplicationType GenerateAtomHybridizationDescriptors_Instance;

      //! @brief Clone function
      //! @return pointer to new GenerateAtomHybridizationDescriptors
      GenerateAtomHybridizationDescriptors *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief initializes the command object for that executable
      util::ShPtr< command::Command> InitializeCommand() const;

      //! @brief returns readme information
      //! @return string containing information about application
      const std::string &GetReadMe() const;

      //! @brief the Main function
      //! @return error code - 0 for success
      int Main() const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM)
      {
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        return OSTREAM;
      }

    //////////////////////
    // helper functions //
    //////////////////////

    private:

      //! @brief write the initializers for chemistry::BondLengths
      //! @param OUTPUT output stream to write to
      void WriteBondLengthsInitializer( std::ostream &OUTPUT) const;

      //! @brief generate the descriptors and bond length info as requested
      //! @param FRAG the fragment for which to calculate the descriptors and bond lengths info
      void Process( const chemistry::FragmentComplete &FRAG) const;

      //! @brief generate counts of element-element bonds for statistics
      //! @param FRAG the fragment for which to calculate the counts
      void CountElementBonds( const chemistry::FragmentEnsemble &FRAG) const;

    }; // GenerateAtomHybridizationDescriptors

  } // namespace app
} // namespace bcl

#endif // BCL_APP_GENERATE_ATOM_HYBRIDIZATION_DESCRIPTORS_H_

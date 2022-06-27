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

#ifndef BCL_APP_SET_SAMPLE_BY_PARTS_ATOMS_H_
#define BCL_APP_SET_SAMPLE_BY_PARTS_ATOMS_H_

// include the interface for all apps
#include "app/bcl_app_apps.h"
#include "app/bcl_app_interface_release.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_configurational_bond_type_data.h"
#include "chemistry/bcl_chemistry_conformation_graph_converter.h"
#include "chemistry/bcl_chemistry_fragment_feed.h"
#include "chemistry/bcl_chemistry_fragment_map_conformer.h"
#include "chemistry/bcl_chemistry_fragment_split_interface.h"
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter_check_serializable.h"
#include "io/bcl_io_file.h"
#include "storage/bcl_storage_template_instantiations.h"
#include "util/bcl_util_implementation.h"
namespace bcl
{
  namespace app
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SetSampleByPartsAtoms
    //! @brief Assigns atoms to the SampleByParts MDL property based on a reference structure to be used with SampleConformations
    //!
    //! @author brownbp1
    //! @date October 23, 2020
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SetSampleByPartsAtoms :
      public InterfaceRelease
    {
    public:

        //! A static instance of this class
        static const ApplicationType SetSampleByPartsAtoms_Instance;

    private:

    //////////
    // data //
    //////////

      // Flag data

      //! output filename
      util::ShPtr< command::FlagInterface> m_OutputFilenameFlag;

      //! Use this atom type data for comparison
      util::ShPtr< command::FlagInterface> m_AtomTypeDataFlag;

      //! Use this bond type data for comparison
      util::ShPtr< command::FlagInterface> m_BondTypeDataFlag;

      //! Reference fragment
      util::ShPtr< command::FlagInterface> m_ReferenceFragmentFlag;

      //! Return complement atoms
      util::ShPtr< command::FlagInterface> m_DisableComplementFlag;

      // Non-flag data

      //! output file
      mutable io::OFStream m_OutputFile;

      //! the type of atom comparison to perform
      mutable chemistry::ConformationGraphConverter::AtomComparisonTypeEnum m_AtomComparison;

      //! the type of bond comparison to perform
      mutable chemistry::ConfigurationalBondTypeData::DataEnum m_BondComparison;

      mutable chemistry::FragmentComplete m_ReferenceFragment;

      //! our workhorse for actually getting the correct atom indices
      chemistry::FragmentMapConformer m_FragmentMapper;

    ///////////////////////////////////
    // construction and destruction //
    ///////////////////////////////////

      //! default constructor
      SetSampleByPartsAtoms();

      //! copy constructor; skips i/o streams
      SetSampleByPartsAtoms( const SetSampleByPartsAtoms &PARENT);

    public:

      //! @brief Clone function
      //! @return pointer to new FoldProtein
      SetSampleByPartsAtoms *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief get a description for the app
      //! @return a brief (no more than 3 line) description for the application
      std::string GetDescription() const;

      //! @brief returns readme information
      //! @return string containing information about application
      const std::string &GetReadMe() const;

      //! @brief initializes the command object for that executable
      util::ShPtr< command::Command> InitializeCommand() const;

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
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

    private:

      //! @brief writes a single molecule out
      //! @param MOLECULE the molecule or fragment to write
      void Write( chemistry::FragmentComplete &MOLECULE) const;

    }; // SetSampleByPartsAtoms

  } // namespace app
} // namespace bcl

#endif // BCL_APP_SET_SAMPLE_BY_PARTS_ATOMS_H_

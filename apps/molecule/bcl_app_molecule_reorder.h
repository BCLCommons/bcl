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

#ifndef BCL_APP_MOLECULE_REORDER_H_
#define BCL_APP_MOLECULE_REORDER_H_

// include header of this class

// includes from bcl - sorted alphabetically
#include "app/bcl_app_apps.h"
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_command.h"
#include "command/bcl_command_flag_dynamic.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter_check_file_existence.h"
#include "command/bcl_command_parameter_check_ranged.h"
#include "command/bcl_command_parameter_check_serializable.h"
#include "descriptor/bcl_descriptor_base.h"
#include "io/bcl_io_file.h"
#include "linal/bcl_linal_vector_operations.h"
#include "math/bcl_math_statistics.h"
#include "util/bcl_util_si_ptr_list.h"
#include "util/bcl_util_stopwatch.h"

namespace bcl
{
  namespace app
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MoleculeReorder
    //! @brief Application for reordering molecules: sorting by property value, reversing, or randomly
    //!
    //! @see @link example_app_molecule_reorder.cpp @endlink
    //! @author mendenjl
    //! @date May 11, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MoleculeReorder :
      public InterfaceRelease
    {
    public:

        // static instance of this class
        static const ApplicationType MoleculeReorder_Instance;

    private:

    //////////
    // data //
    //////////

      //! input sdf files containing molecules to edit.  This is the only non-optional flag
      util::ShPtr< command::FlagInterface> m_InputFilenamesFlag;

      //! randomly permute the molecules in the ensemble
      util::ShPtr< command::FlagInterface> m_RandomizeFlag;

      //! properties to sort by
      util::ShPtr< command::FlagInterface> m_SortByPropertyFlag;

      //! properties to sort by
      util::ShPtr< command::FlagInterface> m_ReverseFlag;

      //! The maximum number of molecules to store per sdf
      //! If the number that would be written out exceeds this limit, then any files written out be given a numeric extension, e.g.
      //! big_ensemble.sdf would become big_ensemble.0.sdf big_ensemble.1.sdf big_ensemble.2.sdf ...
      util::ShPtr< command::FlagInterface> m_MaxMoleculesPerSDFFlag;

      //! Set this flag to sort atoms in the molecule by descending cahn-ingold-prelog priority
      //! Useful for applications that require that the atoms in the molecule are aligned
      util::ShPtr< command::FlagInterface> m_Canonicalize;

      //! flag for reordering atoms
      util::ShPtr< command::FlagInterface> m_Atoms;

      //! File to write molecules into
      util::ShPtr< command::FlagInterface> m_OutputFilenameFlag;

      //! The maximum number of molecules to write out
      util::ShPtr< command::FlagInterface> m_OutputMaxFlag;

      //! How to explicitly reorder the input ensemble
      util::ShPtr< command::FlagInterface> m_OutputOrderFlag;

    ///////////////////////////////////
    // construction and destruction //
    ///////////////////////////////////

      //! default constructor
      MoleculeReorder();

    public:

      //! @brief Clone function
      //! @return pointer to new FoldProtein
      MoleculeReorder *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief initializes the command object for that executable
      util::ShPtr< command::Command> InitializeCommand() const;

      //! @brief check that all the parameter choices were valid
      //! @return true if all the parameter choices were valid
      bool CheckParametersAreAcceptable() const;

      //! @brief returns readme information
      //! @return string containing information about application
      const std::string &GetReadMe() const;

      //! @brief get a description for the app
      //! @return a brief (no more than 3 line) description for the application
      std::string GetDescription() const;

      //! @brief the Main function
      //! @return error code - 0 for success
      int Main() const;

      //! @brief reorder the atoms in the ensemble
      //! @param ENSEMBLE ensemble in which to reorder the atoms of the molecule
      void ReorderAtoms( chemistry::FragmentEnsemble &ENSEMBLE) const;

      //! @brief reorder the atoms in the ensemble
      //! @param ENSEMBLE ensemble in which to reorder the atoms of the molecule
      void CanonicalizeAtoms( chemistry::FragmentEnsemble &ENSEMBLE) const;

      //! @brief write out an ensemble, taking into account m_MaxMoleculesPerSDFFlag
      //! @param ENSEMBLE the ensemble to write
      void WriteEnsemble( const chemistry::FragmentEnsemble &ENSEMBLE) const;

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

      //! @brief load the small molecule ensembles in and merge them into a single ensemble
      //! @return a shared pointer to the combined ensemble
      util::ShPtr< chemistry::FragmentEnsemble> LoadSmallMoleculeEnsembles() const;

      //! @brief function that sorts the ensemble by the numeric value of a particular property
      //! @param ENSEMBLE ensemble to sort
      //! @param PROPERTY property to sort by, should return only 1 value
      void SortByNumericProperty
      (
        chemistry::FragmentEnsemble &ENSEMBLE,
        descriptor::CheminfoProperty PROPERTY
      ) const;

    }; // MoleculeReorder

  } // namespace app
} // namespace bcl

#endif // BCL_APP_MOLECULE_REORDER_H_

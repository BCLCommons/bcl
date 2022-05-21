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

#ifndef BCL_APP_MOLECULE_MUTATE_H_
#define BCL_APP_MOLECULE_MUTATE_H_

// include the interface for all apps
#include "app/bcl_app_apps.h"
#include "app/bcl_app_interface.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_fragment_complete.h"
#include "chemistry/bcl_chemistry_fragment_feed.h"
#include "chemistry/bcl_chemistry_fragment_mutate_interface.h"
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter_check_serializable.h"
#include "io/bcl_io_file.h"
#include "math/bcl_math_mutate_interface.h"
#include "storage/bcl_storage_template_instantiations.h"
#include "util/bcl_util_implementation.h"
namespace bcl
{
  namespace app
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MoleculeMutate
    //! @brief Application for editing small molecules
    //! @details Applies mutations to perturb the chemical structures of small molecules. Utilizes the
    //! FragmentMutateInterface to perform mutations, which returns a single molecule per mutate.
    //!
    //! @see @link example_app_molecule_mutate.cpp @endlink
    //! @author brownbp1
    //! @date January 15, 2021
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MoleculeMutate :
      public InterfaceRelease
    {
    public:

        //! A static instance of this class
        static const ApplicationType MoleculeMutate_Instance;

    private:

    //////////
    // data //
    //////////

      //! file indicating implementation of chemistry::FragmentMutateInterface
      util::ShPtr< command::FlagInterface> m_ImplementationFlag;

      //! accumulate all mutations onto a single molecule
      util::ShPtr< command::FlagInterface> m_AccumulateFlag;

      //! whether to recenter the molecules
      util::ShPtr< command::FlagInterface> m_RecenterFlag;

      //! output filename
      util::ShPtr< command::FlagInterface> m_OutputFilenameFlag;

      //! atom indices that can be mutated
      util::ShPtr< command::FlagInterface> m_MutableAtomsFlag;

      mutable util::Implementation< chemistry::FragmentMutateInterface> m_Mutate; //!< obtains a implementation

      //! molecule index currently being mutated
      mutable size_t m_MoleculeIndex;

      //! output file
      mutable io::OFStream m_OutputFile;

    ///////////////////////////////////
    // construction and destruction //
    ///////////////////////////////////

      //! default constructor
      MoleculeMutate();

      //! copy constructor; skips i/o streams
      MoleculeMutate( const MoleculeMutate &PARENT);

    public:

      //! @brief Clone function
      //! @return pointer to new FoldProtein
      MoleculeMutate *Clone() const;

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

    }; // MoleculeMutate

  } // namespace app
} // namespace bcl

#endif // BCL_APP_MOLECULE_MUTATE_H_

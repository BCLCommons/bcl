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

#ifndef BCL_APP_GENERATE_PARTIAL_CHARGES_FILE_H_
#define BCL_APP_GENERATE_PARTIAL_CHARGES_FILE_H_
// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "app/bcl_app.h"
#include "app/bcl_app_interface.h"
#include "chemistry/bcl_chemistry_fragment_mutate_add_med_chem.h"
#include "command/bcl_command_command.h"
#include "command/bcl_command_flag_interface.h"
#include "linal/bcl_linal_vector.h"
#include "util/bcl_util_sh_ptr.h"
#include "util/bcl_util_static_initialization_fiasco_finder.h"

// include headers from the bcl - sorted alphabetically

namespace bcl
{
  namespace app
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class GeneratePartialChargesFile
    //! @brief Compute atomic partial charges for molecules and format the output for use
    //! with Rosetta's molecule_to_params scripts.
    //!
    //! @author brownbp1
    //! @date 04/23/2023
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class GeneratePartialChargesFile :
      public InterfaceRelease
    {

    private:

    //////////
    // data //
    //////////

      //! flag that specifies the output filename
      util::ShPtr< command::FlagInterface> m_OutputFlag;

      //! flag to specify the partial charge types
      util::ShPtr< command::FlagInterface> m_PartialChargeTypeFlag;

      //! charge types for atomic charge assignment
      enum PartialChargeType
      {
        e_SigmaCharge = 0,        //!< Gasteiger partial charges
        e_PiCharge = 1,           //!< Pi charges
        e_TotalCharge = 2,        //!< Sum of sigma and pi charges
        e_VCharge = 3,            //!< VeraChem partial charges
        e_VCharge2 = 4,           //!< VeraChem partial charges with correction
        s_NumberPartialChargeTypes
      };

      //! partial charge type
      mutable PartialChargeType m_PartialChargeType = e_TotalCharge; // default to TotalCharge for backwards-compatibility

    ///////////////////////////////////
    // construction and destruction //
    ///////////////////////////////////

      //! @brief Default constructor
      GeneratePartialChargesFile();

    public:

      //! @brief Clone function
      //! @return pointer to new GeneratePartialChargesFile
      GeneratePartialChargesFile *Clone() const
      {
        return new GeneratePartialChargesFile( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief initializes the command object for this application
      //! @return a ShPtr to a Command containing all of this applications parameters
      util::ShPtr< command::Command> InitializeCommand() const;

      //! @brief get a description for the app
      //! @return a brief (no more than 3 line) description for the application
      std::string GetDescription() const;

      //! @brief returns readme information
      //! @return string containing information about application
      const std::string &GetReadMe() const;

      //! @brief the Main function
      //! @return 0 for success
      int Main() const;

      //! @brief Compute partial charges of MOL
      //! @param MOL the molecule for which atomic partial charges will be computed
      //! @return atomic partial charges of MOL
      const linal::Vector< float> ComputeAtomicPartialCharges( const chemistry::FragmentComplete &MOL) const;

      //! @brief write the file that contains partial charge and element type for each atom in NCAA
      //! @parameter MOL_INDEX: index of the NCAA in the
      //! @parameter PARTIAL_CHARGES: current index of the C backbone atom
      //! @parameter ATOMS: current index of the N backbone atom
      void WritePartialChargeFile
      (
        const size_t &MOL_INDEX,
        const linal::Vector< float> &PARTIAL_CHARGES,
        const chemistry::AtomVector< chemistry::AtomComplete> &ATOMS
      ) const;

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

      // Static instance of GeneratePartialChargesFile
      static const ApplicationType GeneratePartialChargesFile_Instance;

    }; // class GeneratePartialChargesFile

  } // namespace app
} // namespace bcl

#endif // BCL_APP_GENERATE_PARTIAL_CHARGES_FILE_H_

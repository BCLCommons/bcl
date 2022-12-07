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

#ifndef BCL_APP_GENERATE_ROSETTA_NCAA_INSTRUCTIONS_H_
#define BCL_APP_GENERATE_ROSETTA_NCAA_INSTRUCTIONS_H_
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
    //! @class GenerateRosettaNCAAInstructions
    //! @brief Prepare instructions file for generating non-canonical amino acids in Rosetta
    //!
    //! @author brownbp1, vuot2
    //! @date 09/09/2020
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class GenerateRosettaNCAAInstructions :
      public InterfaceRelease
    {

    private:

    //////////
    // data //
    //////////

      //! to append functional groups to backbones
      mutable chemistry::FragmentMutateAddMedChem m_AddMedChem;

      //! flag that specifies the output filename
      util::ShPtr< command::FlagInterface> m_OutputPrefixFlag;

      //! flag that specifies whether
      util::ShPtr< command::FlagInterface> m_Generate3DFlag;

      //! flag that sets sidechain sample by parts
      util::ShPtr< command::FlagInterface> m_SideChainSampleBypartsFlag;

      //! flag to add extra properties
      util::ShPtr< command::FlagInterface> m_ExtraPropertiesFlag;

      //! flag to specify chirality
      util::ShPtr< command::FlagInterface> m_ChiralityFlag;

      //! flag to specify output partial charge files
      util::ShPtr< command::FlagInterface> m_GeneratePartialChargeFileFlag;

      //! flag to specify indices of CA and Chi 1 atoms
      util::ShPtr< command::FlagInterface> m_CaAndChi1IndicesFlag;

    ///////////////////////////////////
    // construction and destruction //
    ///////////////////////////////////

      //! @brief Default constructor
      GenerateRosettaNCAAInstructions();

    public:

      //! @brief Clone function
      //! @return pointer to new GenerateRosettaNCAAInstructions
      GenerateRosettaNCAAInstructions *Clone() const
      {
        return new GenerateRosettaNCAAInstructions( *this);
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

      //! @brief load neutral glycine residue from library
      //! @return the neutral glycine as the ncaa base
      const storage::Pair< bool, chemistry::FragmentComplete> ReadNCAABase() const;

      //! @brief load neutral glycine dipeptide from library as backbone for alpha AA
      //! @return the neutral glycine dipeptide as the ncaa base
      const storage::Triplet< bool, size_t, chemistry::FragmentComplete> ReadDipeptideBackbone
      (
        const std::string &BACKBONE_TYPE
      ) const;

      //! @brief return the chi1 atom indices
      //! @param NCAA: the atom vector of NCAA
      //! @param C_INDEX: the index of backbone C atom
      //! @param N_INDEX: the index of the backbone N atom
      //! @return the chi1 atom index
      const size_t FindChi1Index
      (
        const chemistry::AtomVector< chemistry::AtomComplete> &NCAA,
        const size_t &CA_INDEX,
        const size_t &C_INDEX,
        const size_t &N_INDEX
      ) const;

      //! @brief return the correct type of NCAA backbone (so this can be added to the instruction file) later
      //! @param NCAA: the atom vector of NCAA
      //! @param C_INDEX: the index of backbone C atom
      //! @param N_INDEX: the index of the backbone N atom
      //! @param CHI1_INDEX: the index of the chi1 angle atom
      const std::string FindCAChirarity
      (
        const chemistry::AtomVector< chemistry::AtomComplete> &NCAA,
        const size_t &CA_INDEX,
        const size_t &C_INDEX,
        const size_t &N_INDEX,
        const size_t &CHI1_INDEX
      ) const;

      //! brief write the final Rosetta instructions file
      const std::string WriteRosettaInstructions
      (
        const size_t &NTER_INDEX,
        const size_t &CA_INDEX,
        const size_t &C_INDEX,
        const size_t &O_INDEX,
        const size_t &CHI1_INDEX,
        const storage::Vector< size_t> &IGNORE_INDICES,
        const size_t &UPPER_N_INDEX,
        const size_t &LOWER_C_INDEX,
        const float &FORMAL_CHARGE,
        const std::string &PROPERTIES
      ) const;

      //! brief output the final property list for the SIDECHAIN of NCAA
      //! brief output the final property list for the SIDECHAIN of NCAA
      const std::string GetSidechainPropertiesList
      (
        const float &FORMAL_CHARGE,
        const size_t &RING_NUM,
        const size_t &AROMATIC_NUM,
        const std::string &CA_CHIRARITY,
        const storage::Vector< std::string> &EXTRAS
      ) const;

      //! brief write the file that contains partial charge and element type for each atom in NCAA
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

      // Static instance of GenerateRosettaNCAAInstructions
      static const ApplicationType GenerateRosettaNCAAInstructions_Instance;

    }; // class GenerateRosettaNCAAInstructions

  } // namespace app
} // namespace bcl

#endif // BCL_APP_GENERATE_ROSETTA_NCAA_INSTRUCTIONS_H_

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

#ifndef BCL_FOLD_SETUP_H_
#define BCL_FOLD_SETUP_H_

// include the namespace header
#include "bcl_fold.h"

// include other forward headers - sorted alphabetically
#include "mc/bcl_mc.fwd.hh"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "assemble/bcl_assemble_protein_storage_file.h"
#include "command/bcl_command_flag_interface.h"
#include "opti/bcl_opti_step_status.h"
#include "quality/bcl_quality_measures.h"
#include "quality/bcl_quality_superimpose_measures.h"
#include "storage/bcl_storage_set.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {

    //! @brief flag for returning the fold setup used
    //! @return the fold setup used
    BCL_API const Setup &GetSetup();

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Setup
    //! @brief This class forms as the setup class for any fold protocol related data and flags
    //! @details class has all the flags and data necessary for initializing mutates, scores and various other
    //! data
    //!
    //! @remarks example unnecessary
    //! @author karakam
    //! @date Mar 27, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API Setup :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      //! protein model built from aa sequences to be folded (this is a dummy start model, it is more like a model used in the first step)
      mutable util::ShPtr< assemble::ProteinModel> m_EmptyModel;

      //! native model
      mutable util::ShPtr< assemble::ProteinModel> m_NativeModel;

      //! start model passed by the user
      mutable util::ShPtr< assemble::ProteinModel> m_StartModel;

      //! number of round
      size_t m_NumberRounds;

      //! quality measures
      storage::Set< quality::Measure> m_QualityMeasures;

      //! superimpose measure
      quality::SuperimposeMeasure m_SuperimposeMeasure;

      //! prefix
      std::string m_Prefix;

      //! step statuses
      storage::Set< opti::StepStatusEnum> m_StepStatuses;

      //! protein storage
      util::ShPtr< assemble::ProteinStorageFile> m_Storage;

    public:

      //! @brief construct and return a static instance of the class
      //! @return a static instance of the class
      static const Setup &GetStaticInstance();

      //! @brief reset the static instance
      static void InitializeStaticInstance();

      //! @brief return vector of all flags applicable
      //! @return vector of all flags applicable
      static util::ShPtrVector< command::FlagInterface> GetAllFlags();

      //! @brief return string containing readme information
      //! @return string containing readme information
      static const std::string &GetReadme();

    private:

      //! @brief construct and return a non-const static instance of the class
      //! @return a non-const static instance of the class
      static Setup &GetStaticInstanceNonConst();

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief private default constructor
      Setup();

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief Clone function
      //! @return pointer to new Setup
      Setup *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief return the protein model constructed from aa sequences (given as fastas) that are supposed to be folded
      //! @return the protein model constructed from aa sequences (given as fastas) that are supposed to be folded
      const util::ShPtr< assemble::ProteinModel> &GetEmptyModel() const;

      //! @brief returns the native/template model if one was provided
      //! @return the native/template model if one was provided
      const util::ShPtr< assemble::ProteinModel> &GetNativeModel() const;

      //! @brief return the starting model after adding sses from -start_model file to it
      //!        (this is used to refine an existing starting model)
      //! @return return the starting model after adding sses from -start_model file to it
      const util::ShPtr< assemble::ProteinModel> GetStartModel() const;

      //! @brief number of rounds
      //! @return number of rounds
      size_t GetNumberRounds() const
      {
        return m_NumberRounds;
      }

      //! @brief return quality measures to be calculated
      //! @return quality measures to be calculated
      const storage::Set< quality::Measure> &GetQualityMeasures() const;

      //! @brief return superimpose measure
      //! @return quality superimpose measure
      const quality::SuperimposeMeasure &GetSuperimposeMeasure() const;

      //! @brief gives the prefix object
      //! @return the prefix that is prepended to output files
      const std::string &GetPrefix() const;

      //! @brief gives the set of step statuses that will be printed by printers
      //! @return set of step statuses that will be printed by printers
      const storage::Set< opti::StepStatusEnum> &GetStepStatuses() const;

      //! @brief get the protein storage
      //! @return protein storage
      const util::ShPtr< assemble::ProteinStorageFile> &GetStorage() const;

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

    }; // class Setup

  } // namespace fold
} // namespace bcl

#endif // BCL_FOLD_SETUP_H_

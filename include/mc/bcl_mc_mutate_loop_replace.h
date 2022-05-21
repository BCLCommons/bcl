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

#ifndef BCL_MC_MUTATE_LOOP_REPLACE_H_
#define BCL_MC_MUTATE_LOOP_REPLACE_H_

// include the namespace header
#include "bcl_mc.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "fold/bcl_fold_locator_loop.h"
#include "fold/bcl_fold_loop_library.h"
#include "math/bcl_math_mutate_result.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace mc
  {

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MutateLoopReplace
    //! @brief Adds a loop to the given protein model
    //! @details Uses conformational hashing to add a loop to the given protein model.
    //!
    //! @see @link example_mc_mutate_loop_replace.cpp @endlink
    //! @author fischea
    //! @date Mar 8, 2016
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    class BCL_API MutateLoopReplace :
      public math::MutateInterface< assemble::ProteinModel>
    {

    //////////
    // data //
    //////////

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    private:

      //! locator for unbuilt loops in a given protein model
      fold::LocatorLoop m_LoopLocator;

      //! path to the loop template library
      std::string m_LoopLibraryFilename;

      //! loop template library
      util::ShPtr< fold::LoopLibrary> m_LoopLibrary;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

    public:

      //! @brief construct from members
      //! @param LOOP_LIBRARY_FILENAME path of the loop template library
      MutateLoopReplace( const std::string &LOOP_LIBRARY_FILENAME = GetDefaultLibraryPath());

      //! @brief copy constructor
      //! @return pointer to a new MutateLoopReplace
      MutateLoopReplace *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns the name of this class
      //! @return the name of this class
      const std::string &GetClassIdentifier() const;

      //! @brief returns the default path of the loop template library
      //! @return the default path of the loop template library
      static const std::string &GetDefaultLibraryPath();

      //! @brief get a short name for this class
      //! @return a short name for this class
      const std::string &GetAlias() const;

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief Set the members of this object from the given LABEL
      //! @param LABEL the label containing members that should be read of this class
      //! @param ERROR_STREAM stream with which to write errors
      //! @return return code indicating success or failure
      bool ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM);

    ///////////////
    // operators //
    ///////////////

      //! @brief applies a mutation to the given protein model
      //! @param MODEL protein model to which to apply the mutation
      //! @return result mutating the given protein model
      math::MutateResult< assemble::ProteinModel> operator()( const assemble::ProteinModel &MODEL) const;

    //////////////////////
    // input and output //
    //////////////////////

    }; // class MutateLoopReplace

  } // namespace mc
} // namespace bcl

#endif // BCL_MC_MUTATE_LOOP_REPLACE_H_

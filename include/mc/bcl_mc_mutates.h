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

#ifndef BCL_MC_MUTATES_H_
#define BCL_MC_MUTATES_H_

// include the namespace header
#include "bcl_mc.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "fold/bcl_fold_default_mutates.h"
#include "util/bcl_util_object_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace mc
  {

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Mutates
    //! @brief Initializes all mutates used for protein structure prediction and adds them to the enumerate
    //!
    //! @see @link example_mc_mutates.cpp @endlink
    //! @author fischea
    //! @date Aug 17, 2016
    //!
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////
    class BCL_API Mutates :
      public util::ObjectInterface
    {

    //////////
    // data //
    //////////

    public:

      //! single instance of this class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////
    // loop hashing //
    //////////////////

      //! mutate for adding loops
      fold::Mutate e_MutateLoopHashAdd;

      //! mutate for replacing loops
      fold::Mutate e_MutateLoopHashReplace;

      //! mutate for replacing loops
      fold::Mutate e_MutateLoopHashResize;

      //! mutate for replacing loops
      fold::Mutate e_MutateLoopHashRemove;

    ////////////////////////////////////
    //  construction and destruction  //
    ////////////////////////////////////

    private:

      //! @brief default constructor
      Mutates();

    public:

      //! @brief clone function
      //! @return pointer to a new Mutates
      Mutates *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @returns single instance of this class
      //! @return single instance of this class
      static Mutates &GetInstance();

      //! @brief returns the name of this class
      //! @return the name of this class
      const std::string &GetClassIdentifier() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief initializes the mutates and adds them to the enumerator
      void Initialize();

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief reads members from a given input stream
      //! @param ISTREAM input stream to read members from
      //! @return input stream which members were read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief writes members into a given output stream
      //! @param OSTREAM output stream to write members into
      //! @param INDENT number of indentations
      //! @return output stream into which members were written
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

    private:

      //! @brief initialize mutates required for loop hashing
      void InitializeLoopHashMutates();

    }; // class Mutates

  } // namespace mc
} // namespace bcl

#endif // BCL_MC_MUTATES_H_

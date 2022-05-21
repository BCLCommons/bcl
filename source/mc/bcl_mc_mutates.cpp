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

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "mc/bcl_mc_mutates.h"

// includes from bcl - sorted alphabetically
#include "mc/bcl_mc_mutate_loop_add.h"
#include "mc/bcl_mc_mutate_loop_add_resize.h"
#include "mc/bcl_mc_mutate_loop_remove.h"
#include "mc/bcl_mc_mutate_loop_replace.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace mc
  {

  //////////
  // data //
  //////////

    //! single instance of this class
    const util::SiPtr< const util::ObjectInterface> Mutates::s_Instance
    (
      GetObjectInstances().AddInstance( new Mutates())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    Mutates::Mutates()
    {
    }

    //! @brief copy constructor
    //! @return pointer to a new Mutates
    Mutates *Mutates::Clone() const
    {
      return new Mutates( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @returns single instance of this class
    //! @return single instance of this class
    Mutates &Mutates::GetInstance()
    {
      static Mutates s_single_instance;
      return s_single_instance;
    }

    //! @brief returns the name of this class
    //! @return the name of this class
    const std::string &Mutates::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief initializes the mutates and adds them to the enumerator
    void Mutates::Initialize()
    {
      // initialize the default mutates
      fold::DefaultMutates::GetInstance().InitializeMutates();

      // initialize mutates for loop hashing
      InitializeLoopHashMutates();
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief reads members from a given input stream
    //! @param ISTREAM input stream to read members from
    //! @return input stream which members were read from
    std::istream &Mutates::Read( std::istream &ISTREAM)
    {
      // return the stream
      return ISTREAM;
    }

    //! @brief writes members into a given output stream
    //! @param OSTREAM output stream to write members into
    //! @param INDENT number of indentations
    //! @return output stream into which members were written
    std::ostream &Mutates::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief initialize mutates required for loop hashing
    void Mutates::InitializeLoopHashMutates()
    {
      // get the enumerated mutates
      fold::Mutates &mutates( fold::GetMutates());

      // add a mutate for adding loops
      e_MutateLoopHashAdd = mutates.AddMutate( MutateLoopAdd());

      // add a mutate for replacing loops
      e_MutateLoopHashReplace = mutates.AddMutate( MutateLoopReplace());

      // add a mutate for resizing loops
      e_MutateLoopHashResize = mutates.AddMutate( MutateLoopAddResize());

      // add a mutate for removing loops
      e_MutateLoopHashRemove = mutates.AddMutate( MutateLoopRemove());
    }

  } // namespace mc
} // namespace bcl

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

// include example header
#include "example.h"
// include the header of the class which this example is for
#include "chemistry/bcl_chemistry_reaction_ensemble.h"

// includes from bcl - sorted alphabetically
#include "descriptor/bcl_descriptor_cheminfo_properties.h"
#include "io/bcl_io_file.h"
#include "sdf/bcl_sdf_fragment_factory.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_chemistry_reaction_ensemble.cpp
  //!
  //! @author geanesar
  //! @date Jul 30, 2015
  //! @remarks status complete
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleChemistryReactionEnsemble :
    public ExampleInterface
  {
  public:

    ExampleChemistryReactionEnsemble *Clone() const
    {
      return new ExampleChemistryReactionEnsemble( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    int Run() const
    {
      // load an ensemble directly from an input stream
      io::IFStream input;
      BCL_ExampleMustOpenInputFile( input, AddExampleInputPathToFilename( e_Chemistry, "multi_reaction_1.rxn"));
      chemistry::ReactionEnsemble rxn_ensemble( input);
      io::File::CloseClearFStream( input);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // copy constructor
      util::ShPtr< chemistry::ReactionEnsemble> sp_ensemble( new chemistry::ReactionEnsemble( rxn_ensemble));
      BCL_ExampleCheck( sp_ensemble.IsDefined(), true);
      BCL_ExampleCheck( sp_ensemble->GetSize() == rxn_ensemble.GetSize(), true);

    /////////////////
    // data access //
    /////////////////

      // check static class name
      BCL_ExampleCheck
      (
        GetStaticClassName< chemistry::ReactionEnsemble>(),
        rxn_ensemble.GetClassIdentifier()
      );

      BCL_ExampleCheck( rxn_ensemble.IsEmpty(), false);
      BCL_ExampleCheck( rxn_ensemble.GetSize(), 1);

      // Check that iterators work
      chemistry::ReactionEnsemble new_ensemble;
      for
      (
        chemistry::ReactionEnsemble::const_iterator itr_rxn( rxn_ensemble.Begin()), itr_rxn_end( rxn_ensemble.End());
        itr_rxn != itr_rxn_end;
        ++itr_rxn
      )
      {
        new_ensemble.PushBack( *itr_rxn);
      }

      BCL_ExampleIndirectCheck( new_ensemble.GetSize(), rxn_ensemble.GetSize(), "iterators work");

      size_t rxn_ens_size( rxn_ensemble.GetSize());
      rxn_ensemble.Append( new_ensemble);
      BCL_ExampleIndirectCheck( rxn_ensemble.GetSize(), 2 * rxn_ens_size, "appending reactions works");

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

      // test read write function
      std::stringstream ens_ss;
      ens_ss << new_ensemble;
      chemistry::ReactionEnsemble ensemble_read_in;
      ens_ss >> ensemble_read_in;

      // check whether read in ensemble has the right number of molecules
      BCL_ExampleIndirectCheck( ensemble_read_in.GetSize(), new_ensemble.GetSize(), "I/O is consistent/symmetric");

    //////////////////////
    // helper functions //
    //////////////////////
      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleChemistryReactionEnsemble

  const ExampleClass::EnumType ExampleChemistryReactionEnsemble::s_Instance
  (
    GetExamples().AddEnum( ExampleChemistryReactionEnsemble())
  );

} // namespace bcl

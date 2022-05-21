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
#include "score/bcl_score_phi_psi_with_sspred.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_locator_sse.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "biol/bcl_biol_membrane.h"
#include "math/bcl_math_bicubic_spline.h"
#include "sspred/bcl_sspred_method_handler.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_score_phi_psi_with_sspred.cpp
  //!
  //! @author alexanns
  //! @date Jan 06, 2011
  //! @remarks status complete
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleScorePhiPsiWithSSPred :
    public ExampleInterface
  {
  public:

    ExampleScorePhiPsiWithSSPred *Clone() const
    {
      return new ExampleScorePhiPsiWithSSPred( *this);
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
      // initialize pdb filename to be read
      const std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "1IE9_shifted.pdb"));

      // get model
      assemble::ProteinModel native_model( Proteins::GetModel( pdb_filename));

      // make a copy of the ShPtr to the sequence of chain A
      util::ShPtr< biol::AASequence> sp_sequence( native_model.GetChain( 'A')->GetSequence());

      //create set of SS methods to use for score
      // storage::Set< sspred::Method> ss_methods( sspred::GetMethods().e_JUFO);
      storage::Set< sspred::Method> ss_methods( sspred::GetMethods().e_JUFO, sspred::GetMethods().e_PSIPRED);

      sspred::MethodHandler::ReadPredictionsForProteinModel
      (
        ss_methods, native_model, "1IE9", AddExampleInputPathToFilename( e_Biology, "")
      );

      // create the locator that finds helix 48-56
      assemble::LocatorSSE helix_48_56_locator( 'A', 48, 56);

      // now locate the SSEs from the protein model
      util::SiPtr< const assemble::SSE> sp_helix_48_56( helix_48_56_locator.Locate( native_model));

      // make sure the sse was correctly located
      BCL_ExampleCheck( sp_helix_48_56.IsDefined(), true);

      const biol::Membrane membrane;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      score::PhiPsiWithSSPred def_constr;
      BCL_ExampleCheck( def_constr.GetScheme(), score::PhiPsiWithSSPred::GetDefaultScheme());

      // constructor taking parameters
      score::PhiPsiWithSSPred param_constr( ss_methods);
      BCL_ExampleCheck( param_constr.GetScheme(), score::PhiPsiWithSSPred::GetDefaultScheme());
      BCL_ExampleCheck( param_constr.GetHistogramFilename(), score::PhiPsi::GetDefaultHistogramFilename());
      BCL_ExampleCheck( param_constr.GetEnergyMap().GetSize(), 3);
      BCL_ExampleCheck( param_constr.GetSSMethods().GetSize(), 2);

      // clone constructor
      util::ShPtr< score::PhiPsiWithSSPred> clone_constr( param_constr.Clone());
      BCL_ExampleCheck( clone_constr->GetScheme(), score::PhiPsiWithSSPred::GetDefaultScheme());
      BCL_ExampleCheck( clone_constr->GetHistogramFilename(), score::PhiPsi::GetDefaultHistogramFilename());
      BCL_ExampleCheck( clone_constr->GetEnergyMap().GetSize(), 3);
      BCL_ExampleCheck( clone_constr->GetSSMethods().GetSize(), 2);

    /////////////////
    // data access //
    /////////////////

      // check GetClassIdentifier
      BCL_ExampleCheck
      (
        GetStaticClassName< score::PhiPsiWithSSPred>(), clone_constr->GetClassIdentifier()
      );

      // check GetDefaultScheme
      BCL_ExampleCheck( score::PhiPsiWithSSPred::GetDefaultScheme(), "phipsi_sspred");

      // check GetScheme
      BCL_ExampleCheck( param_constr.GetScheme(), score::PhiPsiWithSSPred::GetDefaultScheme());

      // GetHistogramFilename
      BCL_ExampleCheck( clone_constr->GetHistogramFilename(), score::PhiPsi::GetDefaultHistogramFilename());

      // GetEnergyMap
      BCL_ExampleCheck( clone_constr->GetEnergyMap().GetSize(), 3);

      // GetSSMethods
      BCL_ExampleCheck( clone_constr->GetSSMethods().GetSize(), 2);

    ///////////////
    // operators //
    ///////////////

      storage::Pair< double, size_t> result( param_constr( *sp_helix_48_56, membrane));
      BCL_MessageDbg( "number scored entities is " + util::Format()( result.Second()));
      BCL_ExampleCheckWithinTolerance( result.First(), -14.1758, 0.00001);
      BCL_ExampleCheckWithinTolerance( result.Second(), 7, 0.001);

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

      // write the score and then read it back in
      WriteBCLObject( param_constr);
      score::PhiPsiWithSSPred score_read;
      ReadBCLObject( score_read);
      // test reading and writing
      {
        BCL_ExampleCheck( param_constr.GetScheme(), score::PhiPsiWithSSPred::GetDefaultScheme());
        BCL_ExampleCheck( param_constr.GetHistogramFilename(), score::PhiPsi::GetDefaultHistogramFilename());
        BCL_ExampleCheck( param_constr.GetEnergyMap().GetSize(), 3);
        BCL_ExampleCheck( param_constr.GetSSMethods().GetSize(), 2);
      }

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleScorePhiPsiWithSSPred

  const ExampleClass::EnumType ExampleScorePhiPsiWithSSPred::s_Instance
  (
    GetExamples().AddEnum( ExampleScorePhiPsiWithSSPred())
  );

} // namespace bcl

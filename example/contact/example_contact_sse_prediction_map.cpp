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
#include "contact/bcl_contact_sse_prediction_map.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_contact_sse_prediction_map.cpp
  //!
  //! @author karakam
  //! @date
  //! @remarks status empty
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleContactSSEPredictionMap :
    public ExampleInterface
  {
  public:

    ExampleContactSSEPredictionMap *Clone() const
    { return new ExampleContactSSEPredictionMap( *this);}

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
//      // initialize pdb filename
//      const std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "1ubi.pdb"));
//      biol::AASequence seq( *Proteins::GetModel( pdb_filename, biol::GetAAClasses().e_AACaCb).GetSequences()(0));
//
//      storage::Set< sspred::Method> ss_methods;
//      ss_methods.Insert( sspred::GetMethods().e_JUFO);
//      ss_methods.Insert( sspred::GetMethods().e_PSIPRED);
//      ss_methods.Insert( sspred::GetMethods().e_SAM);
//      sspred::MethodHandler::ReadPredictionsForAASequence( ss_methods, seq, "1ubi", AddExampleInputPathToFilename( e_Biology, ""));
//
//      const std::string blast_filename( AddExampleInputPathToFilename( e_Biology, "1ubiA.ascii"));
//      // read blast profile
//      BCL_ExampleMustOpenInputFile( read, blast_filename);
//      biol::BlastProfileHandler::ReadProfileForAASequence( read, seq);
//      io::File::CloseClearFStream( read);
//
//      //////////////////////
//      // PREDICTIONMAP
//      //////////////////////
//
//      // create the PredictionMap
//      util::ShPtr< contact::PredictionMap> prediction_map;
//      BCL_MessageStd( "Creating the prediction map and get results ");
//      prediction_map = util::ShPtr< contact::PredictionMap>( new contact::PredictionMap( seq));
//      BCL_MessageStd( "PredictionMap has been created");
//
//      // create the sse pool
//      BCL_MessageStd( "Creating the SSEPool from predictions");
//      assemble::ProteinModelPool< biol::AACaCb> model_pool( assemble::SSEFactory().BuildPool( seq));
//
//      util::ShPtrVector< assemble::SSE< biol::AACaCb> > sse_pool;
//
//      for( util::ShPtrVector< assemble::ProteinModel< biol::AACaCb> >::const_iterator model_itr( model_pool.GetData().Begin()), model_itr_end( model_pool.GetData().End()); model_itr != model_itr_end; ++model_itr)
//      {
//        sse_pool.Append( ( *model_itr)->GetData());
//      }
//      BCL_MessageStd( "SSE Pool has " + util::Format()( sse_pool.GetSize()) + " sses!");
//
//      // create the ssecontactmap
//      BCL_MessageStd( "Creating the ssecontactmap");
//      prediction::SSEContactMap< biol::AACaCb> sse_contact_map( prediction_map, sse_pool);
//
//      BCL_MessageStd( "Iterating over every sse pair in the pool and getting sse_contact probabilities");
//      // iterate over every sse
//      for( util::ShPtrVector< assemble::SSE< biol::AACaCb> >::const_iterator sse_itr_a( sse_pool.Begin()), sse_itr_end( sse_pool.End()); sse_itr_a != sse_itr_end; ++sse_itr_a)
//      {
//        // iterate over every other sse
//        for( util::ShPtrVector< assemble::SSE< biol::AACaCb> >::const_iterator sse_itr_b( sse_pool.Begin()), sse_itr_end( sse_pool.End()); sse_itr_b != sse_itr_end; ++sse_itr_b)
//        {
//          if( sse_itr_a != sse_itr_b && !biol::DoOverlap( **sse_itr_a, **sse_itr_b))
//          {
//
//            storage::VectorND< 2, util::SiPtr< const assemble::SSE< biol::AACaCb> > > vector_nd( **sse_itr_a, **sse_itr_b);
//            BCL_MessageStd( util::Format()( biol::g_SSTypesDescription[ ( *sse_itr_a)->GetType()]) + "( "
//                      + util::Format()( ( *sse_itr_a)->GetFirstAA()->GetSeqID()) + " to " + util::Format()( ( *sse_itr_a)->GetLastAA()->GetSeqID()) + ")"
//                      + " and "
//                      + util::Format()( biol::g_SSTypesDescription[ ( *sse_itr_b)->GetType()]) + "( "
//                      + util::Format()( ( *sse_itr_b)->GetFirstAA()->GetSeqID()) + " to " + util::Format()( ( *sse_itr_b)->GetLastAA()->GetSeqID()) + ")"
//                      + " --> "
//                      + util::Format()(sse_contact_map.GetPrediction( **sse_itr_a, **sse_itr_b))
//                      + " "
//                      + util::Format()(sse_contact_map.GetPrediction( vector_nd)));
//          }
//        }
//      }

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleSspredictions

  const ExampleClass::EnumType ExampleContactSSEPredictionMap::s_Instance
  (
    GetExamples().AddEnum( ExampleContactSSEPredictionMap())
  );

} // namespace bcl

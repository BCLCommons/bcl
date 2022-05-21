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
#include "util/bcl_util_voxel_grid.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "assemble/bcl_assemble_voxel_grid_aa.h"
#include "biol/bcl_biol_aa_base.h"
#include "biol/bcl_biol_aa_complete.h"
#include "biol/bcl_biol_atom.h"
#include "fold/bcl_fold_default_scores.h"
#include "fold/bcl_fold_scores.h"
#include "score/bcl_score_aa_pair_clash.h"
#include "score/bcl_score_aa_sequence_pair.h"
#include "score/bcl_score_protein_model_score_sum.h"
#include "score/bcl_score_protein_model_sse_pairs.h"
#include "util/bcl_util_stopwatch.h"

// external includes - sorted alphabetically
#include <ctime>

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_score_octree_clash.cpp
  //!
  //! @author mendenjl
  //! @date Jan 1, 2017
  //! @remarks status empty
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleScoreSlicelistBenchmark :
    public ExampleInterface
  {
  public:

      ExampleScoreSlicelistBenchmark *Clone() const
    {
      return new ExampleScoreSlicelistBenchmark( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! returns class name
    //! the class name as const ref std::string
    const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    int Run() const
    {

      BCL_MessageStd( "Initializing Slicelist Benchmark . . .");

      // instantiate pdb filename
      // std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "1ubi.pdb"));
      // std::string pdb_filename_update( AddExampleInputPathToFilename( e_Biology, "1ubi_clashed.pdb"));

//        std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "1IE9.pdb"));
//        std::string pdb_filename_update( AddExampleInputPathToFilename( e_Biology, "1IE9_more_clashed.pdb"));

      //  std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "2y01A.pdb"));
      //  std::string pdb_filename_update( AddExampleInputPathToFilename( e_Biology, "2y01A_clashed.pdb"));

      storage::Map< biol::SSType, size_t> ssetype_min_size;
      ssetype_min_size[ biol::GetSSTypes().HELIX] = 4;
      ssetype_min_size[ biol::GetSSTypes().STRAND] = 4;

      storage::Vector< std::string> filenames
      (
        storage::Vector< std::string>::Create( "1C1D.pdb", "1ubi.pdb", "1ubi_clashed.pdb", "1IE9.pdb", "1IE9_more_clashed.pdb", "2y01A.pdb", "2y01A_clashed.pdb", "4yua.pdb", "2vz8.pdb")
      );
      util::Time second( 1, 0);

      score::ProteinModelSSEPairs score_aapair_clash
      (
        util::CloneToShPtr( score::AASequencePair( score::AAPairClash::GetInstance(), false)),
        false, // no normalization
        score::ProteinModel::e_Sequence,
        "AA pair distance"
      );

      score::ProteinModelSSEPairs protein_model_sse_pairs
      (
        fold::Scores::WrapCacheSSEPairScore( util::CloneToShPtr( score::AAPairClash()), false),
        false, // no normalization
        score::ProteinModel::e_Sequence,
        "AA pair distance"
      );

      for( storage::Vector< std::string>::const_iterator itr_filenames( filenames.Begin()), itr_filenames_end( filenames.End()); itr_filenames != itr_filenames_end; ++itr_filenames)
      {
        //instantiate proteinmodels of chains
        assemble::ProteinModel model_01
        (
          Proteins::GetModel
          (
            AddExampleInputPathToFilename( e_Biology, *itr_filenames),
            biol::GetAAClasses().e_AAComplete,
            ssetype_min_size
          )
        );

        const util::SiPtrVector< const biol::AABase> aas( model_01.GetAminoAcids());

        BCL_MessageStd( *itr_filenames + " Number of AAs: " + util::Format()( aas.GetSize()));

        storage::Vector< biol::AAComplete> aa_list;
        storage::List< storage::Pair< biol::AAComplete, linal::Vector3D> > aa_datapoint_list;

        for( size_t i = 0; i < aas.GetSize(); ++i)
        {
          aa_list.PushBack( *aas( i));
          storage::Pair< biol::AAComplete, linal::Vector3D> aa_pair( *aas( i), aas( i)->GetFirstSidechainAtom().GetCoordinates());
          aa_datapoint_list.Append( aa_pair);
        }
        // initialize clash score, this ensures that its data files are read in before any benchmarking occurs
        const score::AAPairClash &inst( score::AAPairClash::GetInstance());
        inst.GetClassIdentifier();

//        util::ShPtr< assemble::SlicelistManagerAA> slm( new assemble::SlicelistManagerAA(util::ConvertToSiPtrVector(aa_list)));

//        util::ShPtr< assemble::ProteinModelData> sp_data( model_01.GetProteinModelData());
//        assemble::ProteinModelData model1_data( *sp_data);
//        sp_data->Insert( assemble::ProteinModelData::e_SlicelistManager, slm);
//        model_01.SetProteinModelData( sp_data);

        //
        //      const util::SiPtrVector< const biol::AABase> aas2(model_02.GetAminoAcids());
        //      storage::Vector<  biol::AAComplete> aa_list2;
        //      storage::List< storage::Pair< biol::AAComplete, linal::Vector3D >> aa_datapoint_list2;
        //
        //      for(size_t i = 0; i< aas2.GetSize(); ++i ) {
        //        aa_list2.Append(*aas2(i));
        //        storage::Pair< biol::AAComplete, linal::Vector3D > aa_pair(*aas2(i), aas2(i)->GetFirstSidechainAtom().GetCoordinates());
        //        aa_datapoint_list2.Append(aa_pair);
        //      }
        //
        //      util::ShPtr< assemble::SlicelistManagerAA< biol::AAComplete>> slm_2( new assemble::SlicelistManagerAA< biol::AAComplete>(util::ConvertToSiPtrVector(aa_list2)));
        //
        //      sp_data = model_02.GetProteinModelData();
        //      assemble::ProteinModelData model2_data( *sp_data);
        //      sp_data->Insert( assemble::ProteinModelData::e_SlicelistManager, slm_2);
        //      model_02.SetProteinModelData( sp_data);

        // ---------------- Slicelist --------------------------------------------
        util::ShPtr< assemble::ProteinModel> model_1_temp( model_01.HardCopy());

        //      util::ShPtr< assemble::ProteinModel> model_2_temp(model_02.HardCopy());

        double actual_score( 0.0);

        //std::clock_t start, loop_start;
        //double duration, list_copy_duration( 0.0);

        //start = std::clock();
        //actual_score = 0.0;
        double itr_per_sec( 0);
//        assemble::SlicelistManagerAA slm_assemble;
//        for(size_t i = 0; i < 5000; ++i) {
//          util::SiPtrVector< const biol::AABase> aas( model_1_temp->GetAminoAcids());
//          slm_assemble.Update( aas);
//          actual_score += slm->GetAAClashScore( aas);
//        }
//
//        duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
//        double list_duration(duration);
//        BCL_MessageStd( " Slicelist Duration: " + util::Format()(list_duration));
//        BCL_MessageStd( " Slicelist Score: " + util::Format()(actual_score));
//
//        start = std::clock();
//        actual_score = 0.0;
//        assemble::Slicelist3dAA sl3d_assemble( 4.0);
//        for(size_t i = 0; i < 5000; ++i) {
//          util::SiPtrVector< const biol::AABase> aas( model_1_temp->GetAminoAcids());
//          sl3d_assemble.SetObjects( aas);
//          actual_score += sl3d_assemble.GetAAClashScore();
//        }
//
//        duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
//        BCL_MessageStd( " Slicelist3d Duration: " + util::Format()(duration));
//        BCL_MessageStd( " Slicelist3d Score: " + util::Format()(actual_score));

//        util::SliceListManagerAA< biol::AAComplete> slm_util;
//        start = std::clock();
//        actual_score = 0.0;
//        list_copy_duration = 0.0;
//        for(size_t i = 0; i < 5000; ++i) {
//          slm_util.Update(util::ConvertToSiPtrVector(aa_list));
//          actual_score += slm_util.GetAAClashScore(model_1_temp->GetAminoAcids());
//        }
//
//        list_duration = duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
//        BCL_MessageStd( " Old Slicelist Duration: " + util::Format()(list_duration));
//        BCL_MessageStd( " Old Slicelist Score: " + util::Format()(actual_score));
//
        auto siptr_aalist( util::ConvertToSiPtrVector( aa_list));
        util::Stopwatch watch( false);
        assemble::VoxelGridAA voxelg( inst.GetDistanceCutoff());
        actual_score = 0.0;
        watch.Reset();
        size_t iterations( 0);
        for( watch.Start(); watch.GetProcessDuration() < second;)
        {
          for( int i( 0); i < 10; ++i, ++iterations)
          {
            voxelg.SetObjects( siptr_aalist);
            actual_score = voxelg.GetAAClashScore( siptr_aalist);
          }
        }
        itr_per_sec = double( iterations) / watch.GetProcessDuration().GetSecondsFractional();
        BCL_MessageStd
        (
          " Voxel Grid neighbor-by-neighbor Iterations/Second: "
          + util::Format()( itr_per_sec)
          + " " + util::Format()( 5000.0 / itr_per_sec) + "s per 5000 iterations"
        );
        BCL_MessageStd( " VoxelGrid Score: " + util::Format()( actual_score));
        BCL_MessageStd( " Voxel grid is using " + util::Format()( voxelg.GetDimension()) + " dimensions to cache positions");
        watch.Stop();
        iterations = 0;

        assemble::VoxelGridAA voxelgc( inst.GetDistanceCutoff(), true);
        iterations = 0;
        for( watch.Start(); watch.GetProcessDuration() < second;)
        {
          for( int i( 0); i < 10; ++i, ++iterations)
          {
            voxelgc.SetObjects( siptr_aalist);
            actual_score = voxelgc.GetAAClashScore( siptr_aalist);
          }
        }
        itr_per_sec = double( iterations) / watch.GetProcessDuration().GetSecondsFractional();
        BCL_MessageStd
        (
          " Voxel Grid cached edge neighbor-by-neighbor Iterations/Second: "
          + util::Format()( itr_per_sec)
          + " " + util::Format()( 5000.0 / itr_per_sec) + "s per 5000 iterations"
        );
        BCL_MessageStd( " VoxelGrid Score: " + util::Format()( actual_score));
        watch.Stop();

        actual_score = 0.0;
        iterations = 0;
        for( watch.Start(); watch.GetProcessDuration() < second;)
        {
          for( int i( 0); i < 10; ++i, ++iterations)
          {
            voxelg.SetObjects( siptr_aalist);
            actual_score = voxelg.GetAAClashScore();
          }
        }
        itr_per_sec = double( iterations) / watch.GetProcessDuration().GetSecondsFractional();
        BCL_MessageStd
        (
          " Voxel Grid global Iterations/Second: "
          + util::Format()( itr_per_sec)
          + " " + util::Format()( 5000.0 / itr_per_sec) + "s per 5000 iterations"
        );
        BCL_MessageStd( " VoxelGrid Score: " + util::Format()( actual_score));
        watch.Stop();

        iterations = 0;
        for( watch.Start(); watch.GetProcessDuration() < second;)
        {
          for( int i( 0); i < 10; ++i, ++iterations)
          {
            voxelgc.SetObjects( siptr_aalist);
            actual_score = voxelgc.GetAAClashScore();
          }
        }
        itr_per_sec = double( iterations) / watch.GetProcessDuration().GetSecondsFractional();
        BCL_MessageStd
        (
          " Voxel Grid cached edge global Iterations/Second: "
          + util::Format()( itr_per_sec)
          + " " + util::Format()( 5000.0 / itr_per_sec) + "s per 5000 iterations"
        );
        BCL_MessageStd( " VoxelGrid cached edge Score: " + util::Format()( actual_score));
        watch.Stop();

        iterations = 0;
        for( watch.Start(); watch.GetProcessDuration() < second;)
        {
          for( int i( 0); i < 10; ++i, ++iterations)
          {
            assemble::VoxelGridAA voxelgl( inst.GetDistanceCutoff(), false);
            voxelgl.SetObjects( siptr_aalist);
            actual_score = voxelgl.GetAAClashScore();
          }
        }
        itr_per_sec = double( iterations) / watch.GetProcessDuration().GetSecondsFractional();
        BCL_MessageStd
        (
          " Voxel Grid global Duration, recreate each time: "
          + util::Format()( itr_per_sec)
          + " " + util::Format()( 5000.0 / itr_per_sec) + "s per 5000 iterations"
        );
        watch.Stop();
        BCL_MessageStd( " VoxelGrid Score: " + util::Format()( actual_score));

        iterations = 0;
        for( watch.Start(); watch.GetProcessDuration() < second;)
        {
          for( int i( 0); i < 10; ++i, ++iterations)
          {
            assemble::VoxelGridAA voxelgcl( inst.GetDistanceCutoff(), true);
            voxelgcl.SetObjects( siptr_aalist);
            actual_score = voxelgcl.GetAAClashScore();
          }
        }
        itr_per_sec = double( iterations) / watch.GetProcessDuration().GetSecondsFractional();
        BCL_MessageStd
        (
          " Voxel Grid cached edge global Duration, recreate each time: "
          + util::Format()( itr_per_sec)
          + " " + util::Format()( 5000.0 / itr_per_sec) + "s per 5000 iterations"
        );
        BCL_MessageStd( " VoxelGrid Score: " + util::Format()( actual_score));
        watch.Stop();

        // ---------------- Comparison with Brute Force ------------------------

        assemble::ProteinModel model_brute_01( *model_1_temp);

        util::ShPtr< assemble::ProteinModel> model_1_brute_temp( model_brute_01.HardCopy());

        iterations = 0;
        for( watch.Start(); watch.GetProcessDuration() < second;) {
          for( int i( 0); i < 10; ++i, ++iterations)
          {
            actual_score = score_aapair_clash( *model_1_temp);
          }
        }
        itr_per_sec = double( iterations) / watch.GetProcessDuration().GetSecondsFractional();
        BCL_MessageStd
        (
          " Brute force iterations / second: "
          + util::Format()( itr_per_sec)
          + " " + util::Format()( 5000.0 / itr_per_sec) + "s per 5000 iterations"
        );
        BCL_MessageStd( " Brute force Score: " + util::Format()( actual_score));
        watch.Stop();

        iterations = 0;
        for( watch.Start(); watch.GetProcessDuration() < second;) {
          for( int i( 0); i < 10; ++i, ++iterations)
          {
            actual_score = protein_model_sse_pairs( *model_1_temp);
          }
        }
        itr_per_sec = double( iterations) / watch.GetProcessDuration().GetSecondsFractional();
        BCL_MessageStd
        (
          " Voxel Grid pairs iterations / second: "
          + util::Format()( itr_per_sec)
          + " " + util::Format()( 5000.0 / itr_per_sec) + "s per 5000 iterations"
        );
        BCL_MessageStd( " Voxel Grid pairs Score: " + util::Format()( actual_score));
        watch.Stop();

//
//        double counter(0);
//        double brute_copy_duration(0.0);
//
//        for(size_t i = 0; i < 5000; ++i) {
//
//         /* const util::SiPtrVector< const biol::AABase> temp_aas(model_1_temp->GetAminoAcids());
//
//          for(size_t i = 0; i< temp_aas.GetSize(); ++i ) {
//            for(size_t j = 0; j< temp_aas.GetSize(); ++j ) {
//              if(clash(*temp_aas(i), *temp_aas(j)) > 0 )
//              {
//                ++counter;
//              }
//            }
//          }  */
//
//
//          counter += e_ScoreAAPairClash->operator()( *model_1_brute_temp);
//
//          loop_start = std::clock();
//  //        util::ShPtr< assemble::ProteinModel> model_temp(model_1_temp->HardCopy());
//  //        model_1_temp = util::ShPtr< assemble::ProteinModel>( model_2_temp->HardCopy());
//  //        model_2_temp = model_temp;
//          brute_copy_duration += ( std::clock() - loop_start ) / (double) CLOCKS_PER_SEC;
//        }
//
//        duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
//
//        BCL_MessageStd( " Brute Force Duration: " + util::Format()(duration));
//        BCL_MessageStd( " Loop Duration: " + util::Format()(brute_copy_duration));
//
//
//        BCL_MessageStd( " Slicelist Duration without Hard Copy surplus: " + util::Format()(list_duration - (list_copy_duration - brute_copy_duration)));
//        BCL_MessageStd( " Brute Force Score: " + util::Format()(counter));
      }
      BCL_MessageStd( " Finished Slicelist Benchmark.");
      // End
      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleScoreSlicelistBenchmark

  const ExampleClass::EnumType ExampleScoreSlicelistBenchmark::s_Instance
  (
    GetExamples().AddEnum( ExampleScoreSlicelistBenchmark())
  );

} // namespace bcl


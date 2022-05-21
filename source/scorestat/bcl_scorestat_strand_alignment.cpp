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
// include header for this class
#include "scorestat/bcl_scorestat_strand_alignment.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_collector_sheet.h"
#include "assemble/bcl_assemble_protein_ensemble.h"
#include "assemble/bcl_assemble_protein_model_data.h"
#include "assemble/bcl_assemble_sse.h"
#include "assemble/bcl_assemble_sse_geometry_packer_all_fragment_pairs.h"
#include "biol/bcl_biol_aa.h"
#include "biol/bcl_biol_aa_back_bone_completer.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_serialization.h"
#include "math/bcl_math_histogram_2d.h"
#include "pdb/bcl_pdb_factory.h"
#include "pdb/bcl_pdb_handler.h"
#include "score/bcl_score_aa_pair_sidechain_interaction.h"
#include "score/bcl_score_loop.h"
#include "storage/bcl_storage_pair.h"
#include "util/bcl_util_message.h"
#include "util/bcl_util_wrapper.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace scorestat
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> StrandAlignment::s_Instance
    (
      util::Enumerated< assemble::AnalyzeProteinEnsembleInterface>::AddInstance( new StrandAlignment())
    );

    //! @brief OutputOption as string
    //! @param OUTPUT_OPTION the OutputOption
    //! @return the string for the OutputOption
    const std::string &StrandAlignment::GetOutputOptionName( const OutputOption &OUTPUT_OPTION)
    {
      static const std::string s_Names[] =
      {
        "Table",
        GetStaticClassName< StrandAlignment::OutputOption>()
      };

      return s_Names[ size_t( OUTPUT_OPTION)];
    }

    //! @brief Output filename as string
    //! @param OutputOption the desired Output Type
    //! @return the string for the output file extension
    const std::string &StrandAlignment::GetOutputFileName( const OutputOption &OUTPUT_OPTION)
    {
      static const std::string s_output_filename_extensions[] =
      {
        "strand_alignment.tbl",
        GetStaticClassName< StrandAlignment::OutputOption>()
      };

      return s_output_filename_extensions[ size_t( OUTPUT_OPTION)];
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    StrandAlignment::StrandAlignment() :
        m_OutputOption( e_Table),
        m_VisualizationOutputPath( "./"),
        m_ChainIds( "")
    {
      // nothing to do
    }

    //! @brief virtual copy constructor
    StrandAlignment *StrandAlignment::Clone() const
    {
      return new StrandAlignment( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as constant reference to std::string
    const std::string &StrandAlignment::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gives the string to append to the the end of a filename to identify this analysis
    //! @return gives the string to append to the the end of a filename to identify this analysis
    const std::string &StrandAlignment::GetOutFilePostfix() const
    {
      return GetOutputFileName( m_OutputOption);
    }

    //! @brief gets the cutoff for distance between Carbonyl-Oxygen and Amide_Nitrogen
    //! @return the cutoff for distance between Carbonyl-Oxygen and Amide_Nitrogen
    const double &StrandAlignment::GetHydrogenBondOHCutoff() const
    {
      return m_HydrogenBondOHCutoff;
    }

    //! @brief gets the path where to write the pdbs to visualize measurements
    //! @return the path where to write the pdbs to visualize measurements
    const std::string &StrandAlignment::GetVisualizationOutputPath() const
    {
      return m_VisualizationOutputPath;
    }

    //! @brief gets chain ids
    //! @return chain ids
    const std::string &StrandAlignment::GetChainIds() const
    {
      return m_ChainIds;
    }

    //! @brief gets the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &StrandAlignment::GetAlias() const
    {
      const static std::string s_name( "StrandAlignment");
      return s_name;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief takes an ensemble and does an analysis resulting in string containing formatted information
    //! @param ENSEMBLE the protein ensemble that will be analyzed
    //! @return string which has the analyzed information about the ensemble
    std::string StrandAlignment::operator()( const assemble::ProteinEnsemble &ENSEMBLE) const
    {
      // make sure the protein ensemble is not empty
      BCL_Assert( ENSEMBLE.GetSize() != 0, "protein ensemble is empty");

      // table for output
      storage::Table< double> strand_pair_table
      (
        storage::TableHeader
        (
          storage::Vector< std::string>::Create
          (
            // topology=1,0 (anti)parallel; on=1,0 this line is a O-N distance and C-O-N angle;
            // oh=0,1 this line is a O-H distance and C-O-H angle, on+oh=1 always;
            "topology", "on", "oh", "seq_id_a", "seq_id_b", "angle_rad", "angle_degree", "cos_angle", "distance"
          )
        )
      );

      // iterate over protein ensemble
      for
      (
        util::ShPtrVector< assemble::ProteinModel>::const_iterator
          protein_model_itr( ENSEMBLE.Begin()), protein_model_itr_end( ENSEMBLE.End());
        protein_model_itr != protein_model_itr_end;
        ++protein_model_itr
      )
      {
        // get current pdb name before all loops start
        const util::ShPtr< util::Wrapper< std::string> >
          &model_name_ptr( ( *protein_model_itr)->GetProteinModelData()->GetData( assemble::ProteinModelData::e_PDBFile));
        const std::string &model_name( model_name_ptr.IsDefined() ? model_name_ptr->GetData() : "");

        // make a copy of current protein model
        assemble::ProteinModel protein_model( *( *protein_model_itr)->Clone());

        // create complete model
        biol::AABackBoneCompleter backbone_completer( true, false, false);
        protein_model = *backbone_completer.CompleteProteinModel( protein_model);

        // get all sheets in the current model
        const assemble::CollectorSheet sheet_collector;
        const util::ShPtrVector< assemble::Domain> sheets( sheet_collector.Collect( protein_model));

        size_t sheet_number( 0);

        // iterate over all sheets
        for
        (
          util::ShPtrVector< assemble::Domain>::const_iterator
            sheet_itr( sheets.Begin()), sheet_itr_end( sheets.End());
          sheet_itr != sheet_itr_end;
          ++sheet_itr, ++sheet_number
        )
        {
          // get the current sheet
          const assemble::Domain &current_sheet( **sheet_itr);

          // use the graph representation of the sheet to find interacting strands
          const assemble::Topology::GraphType &current_sheet_graph( current_sheet.GetTopology()->GetGraph());

          // get a SiPtrVector of all sses in current_sheet
          const util::SiPtrVector< const assemble::SSE> current_sheet_sses( current_sheet.GetSSEs());

          // iterate over all sses in current_sheet
          size_t sse_a_number( 0);
          for
          (
            util::SiPtrVector< const assemble::SSE>::const_iterator
              sse_itr_a( current_sheet_sses.Begin()), sse_itr_end( current_sheet_sses.End());
            sse_itr_a != sse_itr_end;
            ++sse_itr_a, ++sse_a_number
          )
          {
            size_t sse_b_number( sse_a_number + 1);

            // iterate over the rest sses
            for
            (
              util::SiPtrVector< const assemble::SSE>::const_iterator sse_itr_b( sse_itr_a + 1);
                sse_itr_b != sse_itr_end;
              ++sse_itr_b, ++sse_b_number
            )
            {
              // skip strand pairs that are not interacting
              if
              (
                !current_sheet_graph.FindEdge
                (
                  *current_sheet_graph.FindVertex( *sse_itr_a), *current_sheet_graph.FindVertex( *sse_itr_b)
                ).IsDefined()
              )
              {
                continue;
              }

              // get the sse pair and packing information
              const assemble::SSE &strand_a( **sse_itr_a), &strand_b( **sse_itr_b);
              const assemble::SSEGeometryPacking strand_packing( strand_a, strand_b);
              const assemble::SSEGeometryPacking::OrientationEnum packing_orientation( strand_packing.GetOrientation());

              // skip undesired chains
              if( !m_ChainIds.empty() && m_ChainIds.find( strand_a.GetChainID()) == std::string::npos)
              {
                BCL_MessageVrb
                (
                  "Skip strand pair with chain IDs: " + util::Format()( strand_a.GetChainID())
                      + util::Format()( strand_b.GetChainID()) + ", the desired chains are: " + m_ChainIds
                );
                continue;
              }

              // skip unknown orientations
              if
              (
                packing_orientation != assemble::SSEGeometryPacking::e_Parallel
                && packing_orientation != assemble::SSEGeometryPacking::e_AntiParallel
              )
              {
                BCL_MessageVrb
                (
                  "Skip the unknown packing orientation: " + packing_orientation.GetString()
                );
                continue;
              }

              // separately (parallel and anti-parallel) collect atom pairs with potential hydrogen bonds
              std::multimap< int, int> parallel_pairs, antiparallel_pairs;

              // iterate over strand_a
              for
              (
                util::ShPtrVector< biol::AABase>::const_iterator
                  aa_itr_a( strand_a.Begin()), aa_itr_a_end( strand_a.End());
                aa_itr_a != aa_itr_a_end;
                ++aa_itr_a
              )
              {
                // iterate over strand_b
                for
                (
                  util::ShPtrVector< biol::AABase>::const_iterator
                    aa_itr_b( strand_b.Begin()), aa_itr_b_end( strand_b.End());
                  aa_itr_b != aa_itr_b_end;
                  ++aa_itr_b
                )
                {
                  // if strand_a and strand_b are anti-parallel
                  if( packing_orientation == assemble::SSEGeometryPacking::e_AntiParallel)
                  {
                    const int seq_id_a( ( **aa_itr_a).GetSeqID()), seq_id_b( ( **aa_itr_b).GetSeqID());
                    antiparallel_pairs.insert( std::make_pair( seq_id_a, seq_id_b));
                    antiparallel_pairs.insert( std::make_pair( seq_id_b, seq_id_a));
                  }
                  // else if strand_a and strand_b are parallel
                  else if( packing_orientation == assemble::SSEGeometryPacking::e_Parallel)
                  {
                    const int seq_id_a( ( **aa_itr_a).GetSeqID());
                    if( aa_itr_b != strand_b.Begin())
                    {
                      const int seq_id_b_prev( ( **( aa_itr_b - 1)).GetSeqID());
                      parallel_pairs.insert( std::make_pair( seq_id_b_prev, seq_id_a));
                    }
                    if( aa_itr_b != strand_b.End() - 1)
                    {
                      const int seq_id_b_next( ( **( aa_itr_b + 1)).GetSeqID());
                      parallel_pairs.insert( std::make_pair( seq_id_a, seq_id_b_next));
                    }
                  }
                } // end of strand_b
              } // end of stand_a

              // collect the returned data
              storage::List< storage::Vector< double> > data;

              // collect atoms to visualize
              storage::List< biol::Atom> atoms_to_visualize;

              typedef storage::List< storage::Pair< util::ShPtr< biol::AABase>, util::ShPtr< biol::AABase> > > HBList;

              // process anti-parallel strands
              HBList hbonds( SelectShortestHBonds( strand_a, strand_b, antiparallel_pairs));

              for
              (
                HBList::const_iterator hb_itr( hbonds.Begin()), hb_itr_end( hbonds.End());
                  hb_itr != hb_itr_end;
                ++hb_itr
              )
              {
                data.Append
                (
                  ComputeStrandData( hb_itr->First(), hb_itr->Second(), assemble::SSEGeometryPacking::e_AntiParallel)
                );
                atoms_to_visualize.Append( CollectStrandAtoms( hb_itr->First(), hb_itr->Second()));
              }

              // process parallel strands
              hbonds = SelectShortestHBonds( strand_a, strand_b, parallel_pairs);

              for
              (
                HBList::const_iterator hb_itr( hbonds.Begin()), hb_itr_end( hbonds.End());
                  hb_itr != hb_itr_end;
                ++hb_itr
              )
              {
                data.Append
                (
                  ComputeStrandData( hb_itr->First(), hb_itr->Second(), assemble::SSEGeometryPacking::e_Parallel)
                );
                atoms_to_visualize.Append( CollectStrandAtoms( hb_itr->First(), hb_itr->Second()));
              }

              // create the strand string
              const std::string strand_str
              (
                "sheet" + util::Format()( sheet_number) + "_" + packing_orientation.GetString()
                + "_strands" + util::Format()( sse_a_number) + "_" + util::Format()( sse_b_number)
              );

              // insert data into table
              if( m_OutputOption == e_Table)
              {
                for
                (
                  storage::List< storage::Vector< double> >::const_iterator
                    data_itr( data.Begin()), data_itr_end( data.End());
                  data_itr != data_itr_end;
                  ++data_itr
                )
                {
                  strand_pair_table.InsertRow( model_name + "_" + strand_str, *data_itr, true);
                }
              }

              // write pdb files with atoms to visualize
              if( !m_VisualizationOutputPath.empty())
              {
                // split model_name into path and file name
                storage::VectorND< 2, std::string> path_and_filename( io::File::SplitToPathAndFileName( model_name));

                // output_path + filename
                std::string output_path_and_filename
                (
                  m_VisualizationOutputPath + util::GetRuntimeEnvironment().GetPathSeperator()
                  + io::File::RemoveLastExtension( path_and_filename.Second()) + "_" + strand_str + "."
                  + io::File::GetLastExtension( path_and_filename.Second())
                );
                BCL_MessageStd( "Writing visualization file " + output_path_and_filename);

                // get all pdb lines from model
                pdb::Factory factory;
                util::ShPtrList< pdb::Line> visualization_file_lines( factory.WriteCompleteModelToPDBLines( protein_model));

                // add atom lines
                biol::AA amino_acid;
                std::size_t atom_number( visualization_file_lines.GetSize());
                for
                (
                  storage::List< biol::Atom>::const_iterator atom_itr( atoms_to_visualize.Begin()), atom_itr_end( atoms_to_visualize.End());
                    atom_itr != atom_itr_end;
                  ++atom_itr, ++atom_number
                )
                {
                  visualization_file_lines.Append( pdb::Factory::WriteAtomToLine( *atom_itr, amino_acid, 'Z', atom_number));
                }

                // create handler and add lines
                pdb::Handler pdb_handler;
                pdb_handler.AppendLines( visualization_file_lines);

                // write visualization pdb
                io::OFStream pdb_write_stream;
                io::File::MustOpenOFStream( pdb_write_stream, output_path_and_filename);
                pdb_handler.WriteLines( pdb_write_stream);
                io::File::CloseClearFStream( pdb_write_stream);

              } // end of writing pdb files for visualization
            } // end of iterating over the rest sses in current sheet
          } // end of iterating over all sses in current sheet
        } // end of iterating over sheets
      } // end of iterating over protein ensemble

      // write statistics
      std::ostringstream stream;
      if( m_OutputOption == e_Table)
      {
        strand_pair_table.WriteFormatted( stream);
      }
      return stream.str();
    } // end of operator()

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer StrandAlignment::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Computes strand alignment statistics."
      );

      parameters.AddInitializer
      (
        "hydrogen_bond_OH_cutoff",
        "maximum length of a hydrogen bond measured between O-H to be considered for pairing strands",
        io::Serialization::GetAgent( &m_HydrogenBondOHCutoff),
        "4.5"
      );

      parameters.AddInitializer
      (
        "write_strand_alignment_visualization_path",
        "path where to write pdb files that visualize the calculation for the strand alignment data",
        io::Serialization::GetAgent( &m_VisualizationOutputPath),
        "./"
      );

      parameters.AddInitializer
      (
        "chain_ids",
        "string of chain ids to be analyzed, if not specified, take all chains",
        io::Serialization::GetAgent( &m_ChainIds),
        ""
      );

      parameters.AddInitializer
      (
        "output",
        "what type of outputs to provide",
        io::Serialization::GetAgent( &m_OutputOption),
        "Table"
      );

      return parameters;
    }

    //! @brief computes angle between the given three atoms; not using the linal functions b/c of Normalize() call
    //! @param ATOM_A the base atom for calculating vector A->B, A->C
    //! @param ATOM_B the first point
    //! @param ATOM_C the second point
    //! @return the angle in radians
    double StrandAlignment::ComputeAngle( const biol::Atom &ATOM_A, const biol::Atom &ATOM_B, const biol::Atom &ATOM_C) const
    {
      const linal::Vector3D vector_a_b( ( ATOM_B.GetCenter() - ATOM_A.GetCenter()).Normalize());
      const linal::Vector3D vector_a_c( ( ATOM_C.GetCenter() - ATOM_A.GetCenter()).Normalize());
      return linal::ProjAngle( vector_a_b, vector_a_c);
    }

    //! @brief computes the H-bond distance between C-O and H-N, returns undefined if H coordinates are not defined
    //! @brief this function use an oriented H-bond definition: C-O ---> H-N
    //! @param AA_A the first AA ShPtr, the O is taken from this
    //! @param AA_B the second AA ShPtr, the H is taken from this
    //! @return the distance; undefined if a shptr is undefined or if no coordinates for the H are present
    double StrandAlignment::ComputeHBondDistance( const util::ShPtr< biol::AABase> &AA_A, const util::ShPtr< biol::AABase> &AA_B) const
    {
      // if any of AA_A or AA_B is undefined, return nan
      if( !AA_A.IsDefined() || !AA_B.IsDefined())
      {
        return util::GetUndefinedDouble();
      }

      const biol::Atom &a_co( AA_A->GetAtom( biol::GetAtomTypes().O));
      const biol::Atom &b_nh( AA_B->GetAtom( biol::GetAtomTypes().H));

      // if the position of the hydrogen atom is not defined, return nan
      if( !b_nh.AllCoordinatesDefined())
      {
        return util::GetUndefinedDouble();
      }

      return linal::Distance( a_co.GetCenter(), b_nh.GetCenter());
    }

    //! @brief computes a list of values for a single H-bond between the two given AAs
    //! @brief H-bond are defined as oriented: C-O ---> H-N, so the C-O is from AA_A, the H-N is from AA_B
    //! @param AA_A the first AA ShPtr, the C-O is taken from this
    //! @param AA_B the second AA ShPtr, the H-N is taken from this
    //! @param ORIENTATION orientation of the sheet
    //! @return a list of vectors of H-bond specific values identical to the final table columns
    storage::List< storage::Vector< double> > StrandAlignment::ComputeStrandData
    (
      const util::ShPtr< biol::AABase> &AA_A,
      const util::ShPtr< biol::AABase> &AA_B,
      const assemble::SSEGeometryPacking::Orientation &ORIENTATION
    ) const
    {
      BCL_MessageDbg( "AA_A = " + util::Format()( AA_A.IsDefined()) + " AA_B = " + util::Format()( AA_B.IsDefined()));

      // a list of vectors to store strand data
      storage::List< storage::Vector< double> > strand_data;

      // skip undefined hydrogen bonds
      if( !AA_A.IsDefined() || !AA_B.IsDefined())
      {
        BCL_MessageStd( "Skip HBond: one or both of the AA ShPtr is not defined.");
        return strand_data;
      }

      // get interesting atoms C, O, H, N
      const biol::Atom aa_ac( AA_A->GetAtom( biol::GetAtomTypes().C));
      const biol::Atom aa_ao( AA_A->GetAtom( biol::GetAtomTypes().O));
      const biol::Atom aa_bh( AA_B->GetAtom( biol::GetAtomTypes().H));
      const biol::Atom aa_bn( AA_B->GetAtom( biol::GetAtomTypes().N));

      // skip undefined hydrogen atoms
      if( !aa_bh.AllCoordinatesDefined())
      {
        BCL_MessageStd( "Skip HBond: the hydrogen atom has no defined coordinates.");
        return strand_data;
      }

      // get "hydrogen bonding lengths"
      const double distance_ao_bh( linal::Distance( aa_ao.GetCenter(), aa_bh.GetCenter()));
      const double distance_ao_bn( linal::Distance( aa_ao.GetCenter(), aa_bn.GetCenter()));

      // only use O ---> H distance cutoff, more precise than O ---> N
      if( distance_ao_bh > m_HydrogenBondOHCutoff)
      {
        BCL_MessageVrb
        (
          "Skip HBond: O to H distance " + util::Format()( distance_ao_bh) + " > cutoff "
          + util::Format()( m_HydrogenBondOHCutoff)
        );
        return strand_data;
      }

      // get "hydrogen bond angles" and their cosine values
      const double angle_ac_ao_bn( ComputeAngle( aa_ao, aa_ac, aa_bn));
      const double cosine_ac_ao_bn( std::cos( angle_ac_ao_bn));
      BCL_MessageVrb
      (
        "HBond between AA SeqIDs " + util::Format()( AA_A->GetSeqID()) + " and " + util::Format()( AA_B->GetSeqID()) + "\n"
        "C--O--N angle = " + util::Format()( angle_ac_ao_bn) + " in radians or "
        + util::Format()( angle_ac_ao_bn / math::g_Pi * 180) + " in degrees \n"
        + " cosine = " + util::Format()( cosine_ac_ao_bn) + " O--N distance = " + util::Format()( distance_ao_bn) + "\n"
        + " orientation = " + util::Format()( ORIENTATION)
      );
      strand_data.PushBack
      (
        storage::Vector< double>::Create
        (
          ORIENTATION, 1, 0, AA_A->GetSeqID(), AA_B->GetSeqID(),
          angle_ac_ao_bn, angle_ac_ao_bn / math::g_Pi * 180, cosine_ac_ao_bn, distance_ao_bn
        )
      );

      const double angle_ac_ao_bh( ComputeAngle( aa_ao, aa_ac, aa_bh));
      const double cosine_ac_ao_bh( std::cos( angle_ac_ao_bh));

      BCL_MessageVrb
      (
        "HBond between AA SeqID " + util::Format()( AA_A->GetSeqID()) + " and " + util::Format()( AA_B->GetSeqID()) + "\n"
        "C--O--H angle = " + util::Format()( angle_ac_ao_bh) + " in radians or "
        + util::Format()( angle_ac_ao_bh / math::g_Pi * 180) + " in degrees \n"
        + " cosine = " + util::Format()( cosine_ac_ao_bh) + " O--H distance = " + util::Format()( distance_ao_bh) + "\n"
        + " orientation = " + util::Format()( ORIENTATION)
      );
      strand_data.PushBack
      (
        storage::Vector< double>::Create
        (
          ORIENTATION, 0, 1, AA_A->GetSeqID(), AA_B->GetSeqID(),
          angle_ac_ao_bh, angle_ac_ao_bh / math::g_Pi * 180, cosine_ac_ao_bh, distance_ao_bh
        )
      );

      return strand_data;
    } // end of ComputeStrandData

    //! @brief collect all atoms that are part of an H-bond; use an oriented H-bond definition: C-O ---> H-N
    //! @param AA_A the first AA ShPtr, the C-O is taken from this
    //! @param AA_B the second AA ShPtr, the H-N is taken from this
    //! @return the list of atoms of this H-bond
    storage::List< biol::Atom> StrandAlignment::CollectStrandAtoms
    (
      const util::ShPtr< biol::AABase> &AA_A, const util::ShPtr< biol::AABase> &AA_B
    ) const
    {
      BCL_MessageDbg( "AA_A = " + util::Format()( AA_A.IsDefined()) + " AA_B = " + util::Format()( AA_B.IsDefined()));

      // a list of atoms
      storage::List< biol::Atom> strand_atoms;

      // skip undefined amino acid residues
      if( !AA_A.IsDefined() || !AA_B.IsDefined())
      {
        return strand_atoms;
      }

      // get interesting atoms C, O, H, N
      const biol::Atom aa_ac( AA_A->GetAtom( biol::GetAtomTypes().C));
      const biol::Atom aa_ao( AA_A->GetAtom( biol::GetAtomTypes().O));
      const biol::Atom aa_bh( AA_B->GetAtom( biol::GetAtomTypes().H));
      const biol::Atom aa_bn( AA_B->GetAtom( biol::GetAtomTypes().N));

      // skip undefined hydrogen atoms
      if( !aa_bh.AllCoordinatesDefined())
      {
        BCL_MessageVrb( "Skip HBond: the hydrogen atom has no defined coordinates.");
        return strand_atoms;
      }

      // get "hydrogen bonding lengths"
      const double distance_ao_bh( biol::Distance( aa_ao, aa_bh));

      // only use O ---> H distance cutoff, more precise than O ---> N
      if( distance_ao_bh > m_HydrogenBondOHCutoff)
      {
        BCL_MessageVrb
        (
          "Skip HBond: O to H distance " + util::Format()( distance_ao_bh) + " > cutoff "
          + util::Format()( m_HydrogenBondOHCutoff)
        );
        return strand_atoms;
      }

      // push those atoms onto the list and return the list
      strand_atoms.PushBack( aa_ac);
      strand_atoms.PushBack( aa_ao);
      strand_atoms.PushBack( aa_bh);
      strand_atoms.PushBack( aa_bn);

      return strand_atoms;
    } // end of CollectStrandAtoms

    //! @brief selects the shortest possible H-bond (below a cutoff) for each AA from a given set of AA pairs
    //! @param STRAND_A first strand, one AA in each AA pair comes from this strand
    //! @param STRAND_B second strand, one AA in each AA pair comes from this strand
    //! @param POSSIBLE_HBONDS map of possible H-bonds (pairs of AA seqids)
    //! @return a list of pairs of AA ShPtr that are the shortest H-bonds below a cutoff
    storage::List< storage::Pair< util::ShPtr< biol::AABase>, util::ShPtr< biol::AABase> > > StrandAlignment::SelectShortestHBonds
    (
      const assemble::SSE &STRAND_A,
      const assemble::SSE &STRAND_B,
      const std::multimap< int, int> &POSSIBLE_HBONDS
    ) const
    {
      // a list of pairs of AA ShPtrs
      storage::List< storage::Pair< util::ShPtr< biol::AABase>, util::ShPtr< biol::AABase> > > shortest_hbonds;

      // iterate over SeqIDs for first amino acid residues
      for
      (
        std::multimap< int, int>::const_iterator
          hbonds_itr_a( POSSIBLE_HBONDS.begin()), hbonds_itr_end( POSSIBLE_HBONDS.end());
        hbonds_itr_a != hbonds_itr_end;
        hbonds_itr_a = POSSIBLE_HBONDS.upper_bound( hbonds_itr_a->first)
      )
      {
        // shortest distance for current key
        double current_key_shortest_distance( util::GetUndefinedDouble());

        // the ShPtr to AA with shortest distance to current key
        util::ShPtr< biol::AABase> current_key_shortest_distance_aa_sp;

        // get the iterator to the amino acid residue whose SeqID is current key, it could be from either STRAND_A or STRAND_B
        biol::AASequence::const_iterator first_aa_itr_a( STRAND_A.FindAABySeqID( hbonds_itr_a->first));
        biol::AASequence::const_iterator first_aa_itr_b( STRAND_B.FindAABySeqID( hbonds_itr_a->first));

        // this has to be &&, b/c the AA will only be in one of the strands, i.e. one of them always finds the end
        if( first_aa_itr_a == STRAND_A.End() && first_aa_itr_b == STRAND_B.End())
        {
          BCL_MessageStd( "Skipping pair, could not find first AA in SSEs")
          continue;
        }

        // get the ShPtr to the amino acid residue whose SeqID is current key
        util::ShPtr< biol::AABase> first_aa_sp( first_aa_itr_a != STRAND_A.End() ? *first_aa_itr_a : *first_aa_itr_b);

        if( !first_aa_sp.IsDefined())
        {
          BCL_MessageStd( "Skipping pair, first AA ShPtr is not defined")
          continue;
        }

        // iterate over SeqIDs for second amino acid residues
        for
        (
          std::multimap< int, int>::const_iterator
            hbonds_itr_b( hbonds_itr_a), hbonds_itr_b_end( POSSIBLE_HBONDS.upper_bound( hbonds_itr_a->first));
          hbonds_itr_b != hbonds_itr_b_end;
          ++hbonds_itr_b
        )
        {
          // get the iterator to the amino acid residue paired with the amino acid residue whose SeqID is current key
          // it could be from either STRAND_A or STRAND_B
          biol::AASequence::const_iterator second_aa_itr_a( STRAND_A.FindAABySeqID( hbonds_itr_b->second));
          biol::AASequence::const_iterator second_aa_itr_b( STRAND_B.FindAABySeqID( hbonds_itr_b->second));

          if( second_aa_itr_a == STRAND_A.End() && second_aa_itr_b == STRAND_B.End())
          {
            BCL_MessageVrb( "Skipping pair, could not find second AA in SSEs")
            continue;
          }

          // get the ShPtr to the amino acid residue paired with the amino acid residue whose SeqID is current key
          util::ShPtr< biol::AABase> second_aa_sp( second_aa_itr_a != STRAND_A.End() ? *second_aa_itr_a : *second_aa_itr_b);

          if( !second_aa_sp.IsDefined())
          {
            BCL_MessageVrb( "Skipping pair, second AA ShPtr is not defined")
            continue;
          }

          // get hydrogen bonding distance
          const double distance( ComputeHBondDistance( first_aa_sp, second_aa_sp));

          BCL_MessageVrb
          (
            "Distance (" + util::Format()( hbonds_itr_b->first)
            + ", " + util::Format()( hbonds_itr_b->second) + ") = " + util::Format()( distance)
          );

          // get the shortest distance
          if( distance < current_key_shortest_distance || !util::IsDefined( current_key_shortest_distance))
          {
            current_key_shortest_distance = distance;
            current_key_shortest_distance_aa_sp = second_aa_sp;
          }
        }

        if( !current_key_shortest_distance_aa_sp.IsDefined())
        {
          BCL_MessageVrb( "current_key_shortest_distance_aa_sp not defined");
          continue;
        }

        BCL_MessageVrb
        (
          "Shortest Distance (" + util::Format()( hbonds_itr_a->first)
          + ", " + util::Format()( current_key_shortest_distance_aa_sp->GetSeqID())
          + ") = " + util::Format()( current_key_shortest_distance)
        );

        shortest_hbonds.PushBack( storage::Pair< util::ShPtr< biol::AABase>, util::ShPtr< biol::AABase> >( first_aa_sp, current_key_shortest_distance_aa_sp));
      }

      return shortest_hbonds;
    } // end of SelectShortestHBonds

  } // namespace scorestat
} // namespace bcl


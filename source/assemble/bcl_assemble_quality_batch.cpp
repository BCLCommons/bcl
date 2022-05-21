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
#include "assemble/bcl_assemble_quality_batch.h"

// includes from bcl - sorted alphabetically
#include "align/bcl_align_alignment_node.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "assemble/bcl_assemble_quality.h"
#include "contact/bcl_contact_order.h"
#include "contact/bcl_contact_recovery.h"
#include "fold/bcl_fold_default_flags.h"
#include "math/bcl_math_sum_function.h"
#include "score/bcl_score_protein_model_topology.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> QualityBatch::s_Instance
    (
      GetObjectInstances().AddInstance( new QualityBatch())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    QualityBatch::QualityBatch() :
      m_Qualities(),
      m_AtomTypes(),
      m_QualitySuffix( ""),
      m_UseReadableNames( false)
    {
    }

    //! @brief constructor from set of qualities and set of atom types
    //! @param QUALITIES set of qualities to calculate
    //! @param ATOM_TYPES set of atom types to be used
    //! @param QUALITY_SUFFIX suffix to add to end of quality names
    //! @param USE_READABLE_NAMES use readable names (for nice table output)
    QualityBatch::QualityBatch
    (
      const storage::Set< quality::Measure> &QUALITIES,
      const storage::Set< biol::AtomType> &ATOM_TYPES,
      const std::string &QUALITY_SUFFIX,
      const bool USE_READABLE_NAMES
    ) :
      m_Qualities( QUALITIES),
      m_AtomTypes( ATOM_TYPES),
      m_QualitySuffix( QUALITY_SUFFIX),
      m_UseReadableNames( USE_READABLE_NAMES)
    {
      // if the given quality has GDT_TS, add 1A,2A,4A,8A automatically
      if( m_Qualities.Contains( quality::GetMeasures().e_GDT_TS))
      {
        m_Qualities.InsertElement( quality::GetMeasures().e_GDT_1A);
        m_Qualities.InsertElement( quality::GetMeasures().e_GDT_2A);
        m_Qualities.InsertElement( quality::GetMeasures().e_GDT_4A);
        m_Qualities.InsertElement( quality::GetMeasures().e_GDT_8A);
      }
    }

    //! @brief Clone function
    //! @return pointer to new QualityBatch
    QualityBatch *QualityBatch::Clone() const
    {
      return new QualityBatch( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &QualityBatch::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief static function to create the names of columns to be generated from the qualities
    //! @param QUALITIES Set of qualities of interest
    //! @return vector of string that correspond to column names
    storage::Vector< std::string> QualityBatch::ColumnNamesFromQualities
    (
      const storage::Set< quality::Measure> &QUALITIES
    )
    {
      // initialize column names
      storage::Vector< std::string> column_names
      (
        storage::Vector< std::string>::Create( "co", "rco_seq", "rco_aas", "cr6", "cr12", "cr24").Append
        (
          storage::Vector< std::string>::Create( "cr6_mcc", "cr12_mcc", "cr24_mcc", "topology_f").Append
          (
            storage::Vector< std::string>::Create( "topology_tpr", "topology_ppv", "completeness")
          )
        )
      );

      // make a copy of the qualities
      storage::Set< quality::Measure> qualities( QUALITIES);

      // if it has GDT_TS
      if( qualities.Contains( quality::GetMeasures().e_GDT_TS))
      {
        // insert 1A, 2A, 4A and 8A
        qualities.InsertElement( quality::GetMeasures().e_GDT_1A);
        qualities.InsertElement( quality::GetMeasures().e_GDT_2A);
        qualities.InsertElement( quality::GetMeasures().e_GDT_4A);
        qualities.InsertElement( quality::GetMeasures().e_GDT_8A);
      }

      // iterate over measures
      for
      (
        storage::Set< quality::Measure>::const_iterator itr( qualities.Begin()), itr_end( qualities.End());
        itr != itr_end; ++itr
      )
      {
        // insert into column names
        column_names.PushBack( itr->GetName());

        // if this is a RMSD measure
        if
        (
          *itr == quality::GetMeasures().e_RMSD ||
          *itr == quality::GetMeasures().e_RMSD_NoSuperimposition ||
          *itr == quality::GetMeasures().e_RMSD_XYSuperimposition
        )
        {
          // insert also this normalized column name
          column_names.PushBack( itr->GetName() + "100");
        }
        if
        (
          *itr >= quality::GetMeasures().e_GDT_TS && *itr <= quality::GetMeasures().e_GDT_8A
        )
        {
          // also insert the the normalized names for these measures
          column_names.PushBack( itr->GetName() + '_' + "norm");
          column_names.PushBack( itr->GetName() + "_sse_norm");
        }
      }

      // end
      return column_names;
    }

    //! @brief operator calculating the batch of qualities for the given model
    //! @param PROTEIN_MODEL protein model for which the qualities will be calculated
    //! @return list of pairs of quality names and values for the given model
    storage::List< storage::Pair< std::string, double> > QualityBatch::operator()
    (
      const ProteinModel &PROTEIN_MODEL
    ) const
    {
      // initialize boolean to hold whether cache should be used depending on the corresponding command line flag
      static const bool cache( true);

      // construct static object for calculating contact orders
      static const contact::Order co( contact::Order::e_Absolute, "co", cache);
      static const contact::Order rco_seq( contact::Order::e_RelativeSequenceLength, "rco_seq", cache);
      static const contact::Order rco_aas( contact::Order::e_RelativeAAsUsed, "rco_aas", cache);

      // construct static objects for calculating contact recovery
      static const contact::Recovery cr6( 6);
      static const contact::Recovery cr12( 12);
      static const contact::Recovery cr24( 24);
      static const contact::Recovery cr6_mcc( 6, math::ContingencyMatrixMeasures::e_MatthewsCorrelationCoefficient);
      static const contact::Recovery cr12_mcc( 12, math::ContingencyMatrixMeasures::e_MatthewsCorrelationCoefficient);
      static const contact::Recovery cr24_mcc( 24, math::ContingencyMatrixMeasures::e_MatthewsCorrelationCoefficient);

      // construct static object for topology score
      static const score::ProteinModelTopology topology_f( score::ProteinModelTopology::score_f);
      static const score::ProteinModelTopology topology_tpr( score::ProteinModelTopology::score_tpr);
      static const score::ProteinModelTopology topology_ppv( score::ProteinModelTopology::score_ppv);

      // initialize the storage list of pairs
      storage::List< storage::Pair< std::string, double> > results_list;

      // insert the contact order values
      results_list.PushBack( storage::Pair< std::string, double>( "co" + m_QualitySuffix, co( PROTEIN_MODEL)));
      results_list.PushBack( storage::Pair< std::string, double>( "rco_seq" + m_QualitySuffix, rco_seq( PROTEIN_MODEL)));
      results_list.PushBack( storage::Pair< std::string, double>( "rco_aas" + m_QualitySuffix, rco_aas( PROTEIN_MODEL)));

      // get the native model data
      util::ShPtr< ProteinModel> sp_native_model
      (
        PROTEIN_MODEL.GetProteinModelData()->GetData( ProteinModelData::e_NativeFilteredModel)
      );

      // calculate and insert contact recoveries and topology scores, if there is not native model then insert nans
      results_list.PushBack
      (
        storage::Pair< std::string, double>
        (
          "cr6" + m_QualitySuffix,
          sp_native_model.IsDefined() ? cr6( *sp_native_model, PROTEIN_MODEL) : util::GetUndefinedDouble()
        )
      );
      results_list.PushBack
      (
        storage::Pair< std::string, double>
        (
          "cr12" + m_QualitySuffix,
          sp_native_model.IsDefined() ? cr12( *sp_native_model, PROTEIN_MODEL) : util::GetUndefinedDouble()
        )
      );
      results_list.PushBack
      (
        storage::Pair< std::string, double>
        (
          "cr24" + m_QualitySuffix,
          sp_native_model.IsDefined() ? cr24( *sp_native_model, PROTEIN_MODEL) : util::GetUndefinedDouble()
        )
      );
      results_list.PushBack
      (
        storage::Pair< std::string, double>
        (
          "cr6_mcc" + m_QualitySuffix,
          sp_native_model.IsDefined() ? cr6_mcc( *sp_native_model, PROTEIN_MODEL) : util::GetUndefinedDouble()
        )
      );
      results_list.PushBack
      (
        storage::Pair< std::string, double>
        (
          "cr12_mcc" + m_QualitySuffix,
          sp_native_model.IsDefined() ? cr12_mcc( *sp_native_model, PROTEIN_MODEL) : util::GetUndefinedDouble()
        )
      );
      results_list.PushBack
      (
        storage::Pair< std::string, double>
        (
          "cr24_mcc" + m_QualitySuffix,
          sp_native_model.IsDefined() ? cr24_mcc( *sp_native_model, PROTEIN_MODEL) : util::GetUndefinedDouble()
        )
      );
      results_list.PushBack
      (
        storage::Pair< std::string, double>
        (
          "topology_f",
          sp_native_model.IsDefined() ? topology_f( *sp_native_model, PROTEIN_MODEL) : util::GetUndefinedDouble()
        )
      );
      results_list.PushBack
      (
        storage::Pair< std::string, double>
        (
          "topology_tpr",
          sp_native_model.IsDefined() ? topology_tpr( *sp_native_model, PROTEIN_MODEL) : util::GetUndefinedDouble()
        )
      );
      results_list.PushBack
      (
        storage::Pair< std::string, double>
        (
          "topology_ppv",
          sp_native_model.IsDefined() ? topology_ppv( *sp_native_model, PROTEIN_MODEL) : util::GetUndefinedDouble()
        )
      );

      // construct alignments for computing RMSDs
      storage::Map< char, util::ShPtr< align::AlignmentInterface< biol::AABase> > > initial_alignments;
      storage::Map< char, util::ShPtr< align::AlignmentInterface< biol::AABase> > > filtered_alignments;

      // get native model and "SSE-only" model
      util::SiPtr< const ProteinModel> native_model
      (
        PROTEIN_MODEL.GetProteinModelData()->GetData( ProteinModelData::e_NativeModel)
      );
      util::SiPtr< const ProteinModel> filtered_model
      (
        PROTEIN_MODEL.GetProteinModelData()->GetData( ProteinModelData::e_NativeFilteredModel)
      );

      // compute alignments
      if( native_model.IsDefined())
      {
        initial_alignments =
          Quality::CreateAlignmentProteinModels( *native_model, PROTEIN_MODEL);
        filtered_alignments =
          Quality::CreateAlignmentProteinModels( *filtered_model, PROTEIN_MODEL);
      }
      else
      {
        initial_alignments = Quality::CreateAlignmentProteinModelSequencesToSSEs( PROTEIN_MODEL);
        filtered_alignments = Quality::CreateAlignmentFromProteinModelSSEs( PROTEIN_MODEL);
      }

      // remove any assignments in the alignment where the template residues are undefined
      // this can happen if PDB has undefined coordinates. These residues have to be removed because
      // if the evaluated model is a Rosetta model, or a BCL::Fold model generated using predicted pool, the model
      // can actually have coordinates for these residues which are undefined in the template model, which in turn
      // can cause the completeness value to be wrong.
      storage::Map< char, util::ShPtr< align::AlignmentInterface< biol::AABase> > > alignments, sse_alignments;

      // iterate over the initial_alignments
      for
      (
        storage::Map< char, util::ShPtr< align::AlignmentInterface< biol::AABase> > >::const_iterator
          align_itr( initial_alignments.Begin()), align_itr_end( initial_alignments.End());
        align_itr != align_itr_end; ++align_itr
      )
      {
        alignments[ align_itr->first] =
          util::ShPtr< align::AlignmentNode< biol::AABase> >
          (
            new align::AlignmentNode< biol::AABase>
            (
              Quality::RemoveUndefinedTemplateAminoAcidsFromAlignment( *align_itr->second, m_AtomTypes)
            )
          );
      }
      // iterate over the initial_alignments
      for
      (
        storage::Map< char, util::ShPtr< align::AlignmentInterface< biol::AABase> > >::const_iterator
          align_itr( filtered_alignments.Begin()), align_itr_end( filtered_alignments.End());
        align_itr != align_itr_end; ++align_itr
      )
      {
        sse_alignments[ align_itr->first] =
          util::ShPtr< align::AlignmentNode< biol::AABase> >
          (
            new align::AlignmentNode< biol::AABase>
            (
              Quality::RemoveUndefinedTemplateAminoAcidsFromAlignment( *align_itr->second, m_AtomTypes)
            )
          );
      }

      // append the qualities to the list
      results_list.Append( QualitiesFromAlignment( alignments, filtered_alignments));

      // end
      return results_list;
    }

    //! @brief operator calculating the batch of qualities for the two models
    //! @param PROTEIN_MODEL_A protein model for which the qualities will be calculated
    //! @param PROTEIN_MODEL_B protein model for which the qualities will be calculated
    //! @return list of pairs of quality names and values for the given models
    storage::List< storage::Pair< std::string, double> > QualityBatch::operator()
    (
      const ProteinModel &PROTEIN_MODEL_A,
      const ProteinModel &PROTEIN_MODEL_B
    ) const
    {
      // get the alignment for the SSE residues
      const storage::Map< char, util::ShPtr< align::AlignmentInterface< biol::AABase> > > initial_alignments
      (
        Quality::CreateAlignmentProteinModels( PROTEIN_MODEL_A, PROTEIN_MODEL_B)
      );

      // remove any assignments in the alignment where the template residues are undefined
      // this can happen if PDB has undefined coordinates. These residues have to be removed because
      // if the evaluated model is a Rosetta model, or a BCL::Fold model generated using predicted pool, the model
      // can actually have coordinates for these residues which are undefined in the template model, which in turn
      // can cause the completeness value to be wrong.
      storage::Map< char, util::ShPtr< align::AlignmentInterface< biol::AABase> > > alignments;

      // iterate over the initial_alignments
      for
      (
        storage::Map< char, util::ShPtr< align::AlignmentInterface< biol::AABase> > >::const_iterator
          align_itr( initial_alignments.Begin()), align_itr_end( initial_alignments.End());
        align_itr != align_itr_end; ++align_itr
      )
      {
        alignments[ align_itr->first] =
          util::ShPtr< align::AlignmentNode< biol::AABase> >
          (
            new align::AlignmentNode< biol::AABase>
            (
              Quality::RemoveUndefinedTemplateAminoAcidsFromAlignment( *align_itr->second, m_AtomTypes)
            )
          );
      }

      return QualitiesFromAlignment( alignments, alignments);
    }

    //! @brief construct a table that has protein distance measures to be append to a score table
    //! @param PROTEIN_MODEL to be used for measure calculations
    //! @return storage table of measures
    storage::Table< double> QualityBatch::ConstructTable
    (
      const ProteinModel &PROTEIN_MODEL
    ) const
    {
      return ConstructTable( operator()( PROTEIN_MODEL));
    }

    //! @brief construct a table that has protein distance measures to be append to a score table
    //! @param PROTEIN_MODEL_A protein model for which the qualities will be calculated
    //! @param PROTEIN_MODEL_B protein model for which the qualities will be calculated
    //! @return storage table of measures
    storage::Table< double> QualityBatch::ConstructTable
    (
      const ProteinModel &PROTEIN_MODEL_A,
      const ProteinModel &PROTEIN_MODEL_B
    ) const
    {
      return ConstructTable( operator()( PROTEIN_MODEL_A, PROTEIN_MODEL_B));
    }

    //! @brief construct a table that has only the completeness estimate without checking for coordinates
    //! @param PROTEIN_MODEL to be used
    //! @return storage table with completeness estimate
    storage::Table< double> QualityBatch::ContructTableWithCompletenessOnly( const ProteinModel &PROTEIN_MODEL) const
    {
      // get the alignments of sequences
      storage::Map< char, util::ShPtr< align::AlignmentInterface< biol::AABase> > >
        align_seq( Quality::CreateAlignmentFromProteinModelSequences( PROTEIN_MODEL));

      // get the alignments of SSEs
      storage::Map< char, util::ShPtr< align::AlignmentInterface< biol::AABase> > >
        align_sses( Quality::CreateAlignmentFromProteinModelSSEs( PROTEIN_MODEL));

      // calculate AA count in sequences and SSEs
      size_t total_count( 0), sse_count( 0); // needed values to calculate the completeness
      for
      (
        storage::Map< char, util::ShPtr< align::AlignmentInterface< biol::AABase> > >::const_iterator
          align_seq_itr( align_seq.Begin()), align_seq_itr_end( align_seq.End()),
          align_sses_itr( align_sses.Begin()), align_sses_itr_end( align_seq.End());
        align_seq_itr != align_seq_itr_end && align_sses_itr != align_sses_itr_end;
        ++align_seq_itr, ++align_sses_itr
      )
      {
        total_count += align_seq_itr->second->GetSize();
        sse_count += align_sses_itr->second->GetSize();
      }

      const double completeness( ( double)( sse_count) / total_count); // calculate completeness

      // initialize the storage list of pairs and insert the completeness value
      storage::List< storage::Pair< std::string, double> > results;
      results.PushBack( storage::Pair< std::string, double>( "completeness_estimate" + m_QualitySuffix, completeness));

      return ConstructTable( results);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &QualityBatch::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Qualities, ISTREAM);
      io::Serialize::Read( m_AtomTypes, ISTREAM);
      io::Serialize::Read( m_QualitySuffix, ISTREAM);
      io::Serialize::Read( m_UseReadableNames, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &QualityBatch::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Qualities, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_AtomTypes, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_QualitySuffix, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_UseReadableNames, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief calculates a qualities from an alignment
    //! @param ALIGNMENTS the alignments from which the qualities will be calculated
    //! @return list of pairs of string indicating the quality and double indicating the quality value
    storage::List< storage::Pair< std::string, double> > QualityBatch::QualitiesFromAlignment
    (
      const storage::Map< char, util::ShPtr< align::AlignmentInterface< biol::AABase> > > &ALIGNMENTS,
      const storage::Map< char, util::ShPtr< align::AlignmentInterface< biol::AABase> > > &SSE_ALIGNMENTS
    ) const
    {
      // initialize coordinates vector
      storage::VectorND< 2, util::SiPtrVector< const linal::Vector3D> > coordinates_pair;

      // initialize counters
      size_t total_ctr( 0);
      size_t pruned_ctr( 0);
      size_t filtered_ctr( 0);

      // iterate through the map
      for
      (
        storage::Map< char, util::ShPtr< align::AlignmentInterface< biol::AABase> > >::const_iterator
          align_itr( ALIGNMENTS.Begin()), sse_alignments( SSE_ALIGNMENTS.Begin()), align_itr_end( ALIGNMENTS.End());
        align_itr != align_itr_end; ++align_itr
      )
      {
        // get the pruned alignment where gaps and amino acids with undefined coordinates are removed
        const align::AlignmentNode< biol::AABase> pruned_alignment
        (
          Quality::RemoveUndefinedAminoAcidsFromAlignment( *( align_itr->second), m_AtomTypes)
        );

        // update counters
        total_ctr += align_itr->second->GetSize();
        pruned_ctr += pruned_alignment.GetSize();
        filtered_ctr += sse_alignments->second->GetSize();

        // get the common coordinates
        const storage::VectorND< 2, util::SiPtrVector< const linal::Vector3D> > this_coordinates_pair
        (
          Quality::CoordinatesFromAlignment( pruned_alignment, m_AtomTypes)
        );

        // update the total coordinates
        coordinates_pair.First().Append( this_coordinates_pair.First());
        coordinates_pair.Second().Append( this_coordinates_pair.Second());
      }

      // initialize the storage list of pairs
      storage::List< storage::Pair< std::string, double> > results_list;

      // get the ratio of completeness
      const double completeness( double( pruned_ctr) / double( total_ctr));
      const double sse_completeness( double( filtered_ctr) / double( total_ctr));
      const double sse_ratio( filtered_ctr ? completeness / sse_completeness : 0.0);

      // insert the completeness
      results_list.PushBack( storage::Pair< std::string, double>( "completeness" + m_QualitySuffix, completeness));

      // remember whether it has GDT_TS
      const bool has_gdt_ts( m_Qualities.Contains( quality::GetMeasures().e_GDT_TS));

      // iterate over measures
      for
      (
        storage::Set< quality::Measure>::const_iterator
          measure_itr( m_Qualities.Begin()), measure_itr_end( m_Qualities.End());
        measure_itr != measure_itr_end; ++measure_itr
      )
      {
        // if a GDT measure
        if
        (
          *measure_itr >= quality::GetMeasures().e_GDT_TS && *measure_itr <= quality::GetMeasures().e_GDT_8A
        )
        {
          // if it has GDT_TS
          if( has_gdt_ts)
          {
            // if gdt_ts do calculation for all sub GDT measures, otherwise nothing needs to be done
            if( *measure_itr == quality::GetMeasures().e_GDT_TS)
            {
              // calculate individual values and then calculate GDT_TS by averaging to avoid redundant calculations
              // since GDT calculation can be relatively expensive
              const double gdt_1A( ( *quality::GetMeasures().e_GDT_1A)->CalculateMeasure( coordinates_pair.First(), coordinates_pair.Second()));
              const double gdt_2A( ( *quality::GetMeasures().e_GDT_2A)->CalculateMeasure( coordinates_pair.First(), coordinates_pair.Second()));
              const double gdt_4A( ( *quality::GetMeasures().e_GDT_4A)->CalculateMeasure( coordinates_pair.First(), coordinates_pair.Second()));
              const double gdt_8A( ( *quality::GetMeasures().e_GDT_8A)->CalculateMeasure( coordinates_pair.First(), coordinates_pair.Second()));
              const double gdt_ts( ( gdt_1A + gdt_2A + gdt_4A + gdt_8A) / 4.0);

              // insert into the results_vector with the name and the value
              results_list.PushBack
              (
                storage::Pair< std::string, double>( quality::GetMeasures().e_GDT_TS->GetName() + m_QualitySuffix, gdt_ts)
              );
              results_list.PushBack
              (
                storage::Pair< std::string, double>( quality::GetMeasures().e_GDT_TS->GetName() + '_' + "norm" + m_QualitySuffix, gdt_ts * completeness)
              );
              results_list.PushBack
              (
                storage::Pair< std::string, double>( quality::GetMeasures().e_GDT_TS->GetName() + '_' + "sse_norm" + m_QualitySuffix, gdt_ts * sse_ratio)
              );
              results_list.PushBack
              (
                storage::Pair< std::string, double>( quality::GetMeasures().e_GDT_1A->GetName() + m_QualitySuffix, gdt_1A)
              );
              results_list.PushBack
              (
                storage::Pair< std::string, double>( quality::GetMeasures().e_GDT_1A->GetName() + '_' + "norm" + m_QualitySuffix, gdt_1A * completeness)
              );
              results_list.PushBack
              (
                storage::Pair< std::string, double>( quality::GetMeasures().e_GDT_1A->GetName() + '_' + "sse_norm" + m_QualitySuffix, gdt_1A * sse_ratio)
              );
              results_list.PushBack
              (
                storage::Pair< std::string, double>( quality::GetMeasures().e_GDT_2A->GetName() + m_QualitySuffix, gdt_2A)
              );
              results_list.PushBack
              (
                storage::Pair< std::string, double>( quality::GetMeasures().e_GDT_2A->GetName() + '_' + "norm" + m_QualitySuffix, gdt_2A * completeness)
              );
              results_list.PushBack
              (
                storage::Pair< std::string, double>( quality::GetMeasures().e_GDT_2A->GetName() + '_' + "sse_norm" + m_QualitySuffix, gdt_2A * sse_ratio)
              );
              results_list.PushBack
              (
                storage::Pair< std::string, double>( quality::GetMeasures().e_GDT_4A->GetName() + m_QualitySuffix, gdt_4A)
              );
              results_list.PushBack
              (
                storage::Pair< std::string, double>( quality::GetMeasures().e_GDT_4A->GetName() + '_' + "norm" + m_QualitySuffix, gdt_4A * completeness)
              );
              results_list.PushBack
              (
                storage::Pair< std::string, double>( quality::GetMeasures().e_GDT_4A->GetName() + '_' + "sse_norm" + m_QualitySuffix, gdt_4A * sse_ratio)
              );
              results_list.PushBack
              (
                storage::Pair< std::string, double>( quality::GetMeasures().e_GDT_8A->GetName() + m_QualitySuffix, gdt_8A)
              );
              results_list.PushBack
              (
                storage::Pair< std::string, double>( quality::GetMeasures().e_GDT_8A->GetName() + '_' + "norm" + m_QualitySuffix, gdt_8A * completeness)
              );
              results_list.PushBack
              (
                storage::Pair< std::string, double>( quality::GetMeasures().e_GDT_8A->GetName() + '_' + "sse_norm" + m_QualitySuffix, gdt_8A * sse_ratio)
              );
            }
          }
          // if no GDT_TS
          else
          {
            // calculate the value and store
            const double value( ( **measure_itr)->CalculateMeasure( coordinates_pair.First(), coordinates_pair.Second()));

            // insert into the results_vector with the name and the value
            results_list.PushBack( storage::Pair< std::string, double>( measure_itr->GetName() + m_QualitySuffix, value));

            // also insert the the normalized version of this measure
            results_list.PushBack
            (
              storage::Pair< std::string, double>
              (
                measure_itr->GetName() + '_' + "norm" + m_QualitySuffix,
                value * completeness
              )
            ); // also insert the the normalized version of this measure
            results_list.PushBack
            (
              storage::Pair< std::string, double>
              (
                measure_itr->GetName() + '_' + "sse_norm" + m_QualitySuffix,
                value * sse_ratio
              )
            );
          }
        }
        // any other measure
        else
        {
          // calculate the value and store
          const double value( ( **measure_itr)->CalculateMeasure( coordinates_pair.First(), coordinates_pair.Second()));

          // insert into the results_vector with the name and the value
          results_list.PushBack( storage::Pair< std::string, double>( measure_itr->GetName() + m_QualitySuffix, value));

          // if this is a RMSD measure
          if
          (
            *measure_itr == quality::GetMeasures().e_RMSD ||
            *measure_itr == quality::GetMeasures().e_RMSD_NoSuperimposition ||
            *measure_itr == quality::GetMeasures().e_RMSD_XYSuperimposition
          )
          {
            // calculate rmsd_100 value
            const double rmsd_100( Quality::RMSD100( value, pruned_ctr));
            BCL_MessageStd
            (
              "Calculating RMSD100=" + util::Format()( rmsd_100)
              + " from RMSD=" + util::Format()( value) + " and amino_acid_count=" + util::Format()( pruned_ctr)
            );

            // insert also this normalized value into the table
            results_list.PushBack
            (
              storage::Pair< std::string, double>( measure_itr->GetName() + "100" + m_QualitySuffix, rmsd_100)
            );
          }
        }
      }

      return results_list;
    }

    //! @brief construct a table that has protein distance measures to be append to a score table
    //! @param QUALITY_LIST list of calculated quality measures and their values
    //! @return storage table of measures
    storage::Table< double> QualityBatch::ConstructTable
    (
      const storage::List< storage::Pair< std::string, double> > &QUALITY_LIST
    ) const
    {
      // construct storage table from HEADER
      storage::Table< double> measures_table
      (
        math::SumFunctionMixin< score::ProteinModel>::GetValueTableVerticalColumnNames()
      );

      // iterate over the measures vector
      for
      (
        storage::List< storage::Pair< std::string, double> >::const_iterator
          measure_itr( QUALITY_LIST.Begin()), measure_itr_end( QUALITY_LIST.End());
        measure_itr != measure_itr_end; ++measure_itr
      )
      {
        // insert the row for this measure
        measures_table.InsertRow
        (
          m_UseReadableNames ? GetFormattedName( measure_itr->First()) : measure_itr->First(),
          storage::Vector< double>::Create( 1.0, measure_itr->Second(), measure_itr->Second())
        );
      }

      // end
      return measures_table;
    }

    //! @brief get formatted name for better readability
    //! @param ROW_NAME original row name
    //! @return formatted name for better readability
     std::string QualityBatch::GetFormattedName( const std::string &ROW_NAME)
    {
      // initialize static map
      static storage::Map< std::string, std::string> s_name_map;

      if( s_name_map.IsEmpty())
      {
        s_name_map[ "RMSD100"] = storage::Table< double>::s_Indentation + "RMSD100";
        s_name_map[ "cr12"] = storage::Table< double>::s_Indentation + "Contact recovery";
      }

      // check the map and return new string if found
      storage::Map< std::string, std::string>::const_iterator find_itr( s_name_map.Find( ROW_NAME));
      if( find_itr != s_name_map.End())
      {
        return find_itr->second;
      }

      // not in map so return original
      return ROW_NAME;
    }

  } // namespace assemble
} // namespace bcl

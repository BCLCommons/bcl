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
#include "assemble/bcl_assemble_quality.h"

// includes from bcl - sorted alphabetically
#include "align/bcl_align_aligner_dp.h"
#include "align/bcl_align_handler_interface.h"
#include "align/bcl_align_sequence.h"
#include "assemble/bcl_assemble_collector_aa_specified.h"
#include "assemble/bcl_assemble_collector_common_aa.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "biol/bcl_biol_align_by_aa_data.h"
#include "biol/bcl_biol_align_by_pdb_id.h"
#include "command/bcl_command_flag_dynamic.h"
#include "command/bcl_command_parameter.h"
#include "function/bcl_function_binary_sum.h"
#include "io/bcl_io_directory_entry.h"
#include "io/bcl_io_file.h"
#include "score/bcl_score_aa_assignments.h"
#include "score/bcl_score_assignment_with_gap.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> Quality::s_Instance
    (
      GetObjectInstances().AddInstance( new Quality())
    );

    //! @brief gives the flag which allows specifying multiple alignment files corresponding to multiple chains
    //! @return flag which allows specifying multiple alignment files corresponding to multiple chains
    const util::ShPtr< command::FlagInterface> &Quality::GetFlagAlignments()
    {
      static const util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagDynamic
        (
          "alignment",
          "use an alignment for comparing PDB(s) with reference structure",
          command::Parameter
          (
            "alignment_files",
            "list of PIR alignment files, one for each chain, test/template sequence MUST be first - "
            "name must be of format, *<chain_id>.*, "
            "where the chain_id is the char representing which chains are aligned "
            "(i.e. alignA.pir would have the alignment of the A chains)", ""
          )
        )
      );

      return s_flag;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    //! @details construct with rmsd and backbone atom types
    Quality::Quality() :
      m_Measure( quality::GetMeasures().e_RMSD),
      m_AtomTypes( biol::GetAtomTypes().GetBackBoneAtomTypes())
    {
    }

    //! @brief construct from a quality measure and the atom types to be considered
    //! @param QUALITY_MEASURE quality measure to be calculated
    //! @param ATOM_TYPES Set of atom types to be superimposed
    Quality::Quality
    (
      const quality::Measure &QUALITY_MEASURE,
      const storage::Set< biol::AtomType> &ATOM_TYPES
    ) :
      m_Measure( QUALITY_MEASURE),
      m_AtomTypes( ATOM_TYPES)
    {
    }

    //! @brief virtual copy constructor
    //! @return pointer a new Quality copied from this quality
    Quality *Quality::Clone() const
    {
      return new Quality( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Quality::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gets the scheme for this function class
    //! @return the scheme for this function class
    const std::string &Quality::GetScheme() const
    {
      return m_Measure.GetName();
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief calculate the given QUALITY_MEASURE between the coordinates in the SSEs of the Protein model and the original coordinates given from the pdbfile, that are stored in the original sequence
    //! @param QUALITY_MEASURE quality measure to be calculated
    //! @param PROTEIN_MODEL the protein model for which the QUALITY_MEASURE of the coordinates in the SSEs is calculated against the coordinates in the sequence in the chains
    //! @param ATOM_TYPES Set of atom types to be superimposed
    //! @return QUALITY_MEASURE between the coordinates in the SSEs of the Protein model and the original coordinates given from the pdbfile, that are stored in the original sequence
    double Quality::Calculate
    (
      const quality::Measure &QUALITY_MEASURE,
      const ProteinModel &PROTEIN_MODEL,
      const storage::Set< biol::AtomType> &ATOM_TYPES
    )
    {
      // assert that measure is defined
      BCL_Assert( QUALITY_MEASURE.IsDefined(), "Undefined quality measure provided for calculating quality");

      // check that atom types are given
      if( ATOM_TYPES.IsEmpty())
      {
        BCL_MessageCrt( "no atom types given for quality calculation");
        return util::GetUndefined< double>();
      }

      // get the pair of coordinate vectors for each chain that has the current and original coordinates
      storage::Pair< util::SiPtrVector< const linal::Vector3D>, util::SiPtrVector< const linal::Vector3D> >
        model_sequence_coordinates( GetChainAndSequenceCoordinates( PROTEIN_MODEL, ATOM_TYPES));

      // superimposition with less than three coordinates is ambiguous
      if( model_sequence_coordinates.First().GetSize() < 3)
      {
        BCL_MessageCrt( "3 coordinates are needed at least for quality calculation")
        return util::GetUndefined< double>();
      }

      // calculate the quality measure for common coordinates and return it
      return
        ( *QUALITY_MEASURE)->CalculateMeasure
        (
          model_sequence_coordinates.First(), model_sequence_coordinates.Second()
        );
    }

    //! @brief calculate the given QUALITY_MEASURE between the coordinates in the SSEs of the Protein model and the original coordinates given from the pdbfile, that are stored in the original sequence
    //! @param QUALITY_MEASURE quality measure to be calculated
    //! @param PROTEIN_MODEL the protein model for which the QUALITY_MEASURE of the coordinates in the SSEs is calculated against the coordinates in the sequence in the chains
    //! @param NATIVE_MODEL native structure
    //! @param ATOM_TYPES Set of atom types to be superimposed
    //! @return QUALITY_MEASURE between the coordinates in the SSEs of the Protein model and the original coordinates given from the pdbfile, that are stored in the original sequence
    double Quality::Calculate
    (
      const quality::Measure &QUALITY_MEASURE,
      const ProteinModel &PROTEIN_MODEL,
      const ProteinModel &NATIVE_MODEL,
      const storage::Set< biol::AtomType> &ATOM_TYPES
    )
    {
      // assert that measure is defined
      BCL_Assert( QUALITY_MEASURE.IsDefined(), "Undefined quality measure provided for calculating quality");

      // check that atom types are given
      if( ATOM_TYPES.IsEmpty())
      {
        BCL_MessageCrt( "no atom types given for rmsd calculation");
        return util::GetUndefined< double>();
      }

      // collect aas with defined coordinates from protein model
      const util::SiPtrList< const biol::AABase>
        aas_defined_model( CollectorCommonAA::CollectDefinedAAsInSSEs( PROTEIN_MODEL));

      // collect aas with defined coordinates from native model
      const util::SiPtrList< const biol::AABase>
        aas_defined_native( CollectorCommonAA::CollectDefinedAAsInSSEs( NATIVE_MODEL));

      // create VectorND "common_amino_acids"
      storage::VectorND< 2, util::SiPtrList< const biol::AABase> > common_amino_acids
      (
        // initialize with all the amino acids that are common to both "protein_model_a" and "protein_model_b"
        CollectorCommonAA().Collect
        (
          storage::VectorND< 2, util::SiPtr< const util::SiPtrList< const biol::AABase> > >
          (
            util::ToSiPtr( aas_defined_model),
            util::ToSiPtr( aas_defined_native)
          )
        )
      );

      // create VectorND "all_atom_coordinates" to hold the desired atom coordinates of both proteins
      storage::VectorND< 2, util::SiPtrVector< const linal::Vector3D> > all_atom_coordinates
      (
        // initialize with all the atom coordinates
        GetAllDesiredAtomCoordinates( common_amino_acids, ATOM_TYPES)
      );

      // calculate the quality measure for common coordinates and return it
      return ( *QUALITY_MEASURE)->CalculateMeasure( all_atom_coordinates.First(), all_atom_coordinates.Second());
    }

    //! @brief superimpose the coordinates in the SSEs of the Protein model and the original coordinates given from the pdbfile, that are stored in the original sequence
    //! @param SUPERIMPOSE_MEASURE superimpose measure to be used for superimposition
    //! @param PROTEIN_MODEL protein model that will be superimpose to its native coordinates
    //! @param ATOM_TYPES Set of atom types to be superimposed
    //! @return whether superimposition succeeded and the transformation used
    storage::Pair< bool, math::TransformationMatrix3D> Quality::SuperimposeModel
    (
      const quality::SuperimposeMeasure &SUPERIMPOSE_MEASURE,
      ProteinModel &PROTEIN_MODEL,
      const storage::Set< biol::AtomType> &ATOM_TYPES
    )
    {
      // initialize static default transformation
      static const math::TransformationMatrix3D s_default_transformation;

      // assert that measure is defined
      BCL_Assert( SUPERIMPOSE_MEASURE.IsDefined(), "Undefined quality measure provided for calculating superimposition");

      // check that atom types are given
      if( ATOM_TYPES.IsEmpty())
      {
        BCL_MessageCrt( "no atom types given for quality calculation");
        return storage::Pair< bool, math::TransformationMatrix3D>( false, s_default_transformation);
      }

      // get the pair of coordinate vectors for each chain that has the current and original coordinates
      storage::Pair< util::SiPtrVector< const linal::Vector3D>, util::SiPtrVector< const linal::Vector3D> >
        model_sequence_coordinates( GetChainAndSequenceCoordinates( PROTEIN_MODEL, ATOM_TYPES));

      // superimposition with less than three coordinates is ambiguous
      if( model_sequence_coordinates.First().GetSize() < 3)
      {
        BCL_MessageCrt( "3 coordinates are needed at least for quality calculation")
        return storage::Pair< bool, math::TransformationMatrix3D>( false, s_default_transformation);
      }

      // calculate the superimposition transformation matrix
      const math::TransformationMatrix3D transformation
      (
        ( *SUPERIMPOSE_MEASURE)->CalculateSuperimposition
        (
          model_sequence_coordinates.First(), model_sequence_coordinates.Second()
        )
      );

      // if the transformation is not defined
      if( !transformation.IsDefined())
      {
        BCL_MessageCrt( "the superimposition failed");
        return storage::Pair< bool, math::TransformationMatrix3D>( false, s_default_transformation);
      }

      // transform the protein model
      PROTEIN_MODEL.Transform( transformation);

      // return
      return storage::Pair< bool, math::TransformationMatrix3D>( true, transformation);
    }

    //! @brief superimpose given model to another model
    //! @param SUPERIMPOSE_MEASURE superimpose measure to be used for superimposition
    //! @param PROTEIN_MODEL the protein model to be superimposed
    //! @param REFERENCE_MODEL reference structure
    //! @param ATOM_TYPES Set of atom types to be superimposed
    //! @return whether superimposition succeeded and the transformation used
    storage::Pair< bool, math::TransformationMatrix3D> Quality::SuperimposeModel
    (
      const quality::SuperimposeMeasure &SUPERIMPOSE_MEASURE,
      ProteinModel &PROTEIN_MODEL,
      const ProteinModel &REFERENCE_MODEL,
      const storage::Set< biol::AtomType> &ATOM_TYPES
    )
    {
      // initialize static default transformation
      static const math::TransformationMatrix3D s_default_transformation;

      // assert that measure is defined
      BCL_Assert( SUPERIMPOSE_MEASURE.IsDefined(), "Undefined quality measure provided for calculating superimposition");

      // check that atom types are given
      if( ATOM_TYPES.IsEmpty())
      {
        BCL_MessageCrt( "no atom types given for superimposition");
        return storage::Pair< bool, math::TransformationMatrix3D>( false, s_default_transformation);
      }

      // collect aas with defined coordinates from protein models
      const util::SiPtrList< const biol::AABase>
        aas_defined_model( CollectorCommonAA::CollectDefinedAAsInSSEs( PROTEIN_MODEL));
      const util::SiPtrList< const biol::AABase>
        aas_defined_reference( CollectorCommonAA::CollectDefinedAAsInSSEs( REFERENCE_MODEL));

      // create VectorND "common_amino_acids"
      storage::VectorND< 2, util::SiPtrList< const biol::AABase> > common_amino_acids
      (
        // initialize with all the amino acids that are common to both "protein_model" and "reference_model"
        CollectorCommonAA().Collect
        (
          storage::VectorND< 2, util::SiPtr< const util::SiPtrList< const biol::AABase> > >
          (
            util::ToSiPtr( aas_defined_model),
            util::ToSiPtr( aas_defined_reference)
          )
        )
      );

      // create VectorND "all_atom_coordinates" to hold the desired atom coordinates of both proteins
      storage::VectorND< 2, util::SiPtrVector< const linal::Vector3D> > all_atom_coordinates
      (
        // initialize with all the atom coordinates
        GetAllDesiredAtomCoordinates( common_amino_acids, ATOM_TYPES)
      );

      // calculate the superimposition transformation matrix
      const math::TransformationMatrix3D transformation
      (
        ( *SUPERIMPOSE_MEASURE)->CalculateSuperimposition( all_atom_coordinates.First(), all_atom_coordinates.Second())
      );

      // if the transformation is not defined
      if( !transformation.IsDefined())
      {
        BCL_MessageCrt( "the superimposition failed");
        return storage::Pair< bool, math::TransformationMatrix3D>( false, s_default_transformation);
      }

      // transform the protein model
      PROTEIN_MODEL.Transform( transformation);

      // return
      return storage::Pair< bool, math::TransformationMatrix3D>( true, transformation);
    }

    //! @brief superimpose given model to another model
    //! @param SUPERIMPOSE_MEASURE superimpose measure to be used for superimposition
    //! @param PROTEIN_MODEL the protein model to be superimposed
    //! @param REFERENCE_MODEL reference structure
    //! @param ATOM_TYPES Set of atom types to be superimposed
    //! @return whether superimposition succeeded and the transformation used
    storage::Pair< bool, math::TransformationMatrix3D> Quality::SuperimposeModelWithAlignment
    (
      const quality::SuperimposeMeasure &SUPERIMPOSE_MEASURE,
      ProteinModel &PROTEIN_MODEL,
      const ProteinModel &REFERENCE_MODEL,
      const storage::Set< biol::AtomType> &ATOM_TYPES
    )
    {
      // initialize static default transformation
      static const math::TransformationMatrix3D s_default_transformation;

      // assert that measure is defined
      BCL_Assert( SUPERIMPOSE_MEASURE.IsDefined(), "Undefined quality measure provided for calculating superimposition");

      // check that atom types are given
      if( ATOM_TYPES.IsEmpty())
      {
        BCL_MessageCrt( "no atom types given for superimposition");
        return storage::Pair< bool, math::TransformationMatrix3D>( false, s_default_transformation);
      }

      // construct alignments for computing RMSDs
      storage::Map< char, util::ShPtr< align::AlignmentInterface< biol::AABase> > >
        initial_alignments_a, initial_alignments_b, initial_alignments;

      // compute alignments
      initial_alignments_a = Quality::CreateAlignmentFromProteinModelSequences( PROTEIN_MODEL);

      // compute alignments
      initial_alignments_b = Quality::CreateAlignmentFromProteinModelSequences( REFERENCE_MODEL);

      initial_alignments = Quality::AlignSingleAlignments( initial_alignments_a, initial_alignments_b);

      // initialize coordinates vector
      storage::VectorND< 2, util::SiPtrVector< const linal::Vector3D> > coordinates_pair;

      // iterate through the map
      for
      (
        storage::Map< char, util::ShPtr< align::AlignmentInterface< biol::AABase> > >::const_iterator
          align_itr( initial_alignments.Begin()), align_itr_end( initial_alignments.End());
        align_itr != align_itr_end;
        ++align_itr
      )
      {
        // get the pruned alignment where gaps and amino acids with undefined coordinates are removed
        const align::AlignmentNode< biol::AABase> pruned_alignment
        (
          Quality::RemoveUndefinedAminoAcidsFromAlignment( *( align_itr->second), ATOM_TYPES)
        );

        // get the common coordinates
        const storage::VectorND< 2, util::SiPtrVector< const linal::Vector3D> > this_coordinates_pair
        (
          Quality::CoordinatesFromAlignment( pruned_alignment, ATOM_TYPES)
        );

        // update the total coordinates
        coordinates_pair.First().Append( this_coordinates_pair.First());
        coordinates_pair.Second().Append( this_coordinates_pair.Second());
      }

      // calculate the superimposition transformation matrix
      const math::TransformationMatrix3D transformation
      (
        ( *SUPERIMPOSE_MEASURE)->CalculateSuperimposition( coordinates_pair.First(), coordinates_pair.Second())
      );

      // if the transformation is not defined
      if( !transformation.IsDefined())
      {
        BCL_MessageCrt( "the superimposition failed");
        return storage::Pair< bool, math::TransformationMatrix3D>( false, s_default_transformation);
      }

      // transform the protein model
      PROTEIN_MODEL.Transform( transformation);

      // return
      return storage::Pair< bool, math::TransformationMatrix3D>( true, transformation);
    }

    //! @brief calculates RMSD100 by normalizing RMSD by the number of residues
    //! @param RMSD RMSD to be converted to RMSD100
    //! @param NUMBER_RESIDUES number of residues aligned in RMSD calculation
    //! @return RMSD100 calculated from the given RMSD
    double
    Quality::RMSD100
    (
      const double RMSD,
      const size_t NUMBER_RESIDUES
    )
    {
      // if there is less 14 residues (it leads to a negative value)
      if( NUMBER_RESIDUES <= 14)
      {
        BCL_MessageCrt( "There needs to be at least 14 residues for RMSD100 calculation");
        // return undefined
        return util::GetUndefinedDouble();
      }

      // otherwise go ahead with the calculation
      return RMSD / ( 1.0 + std::log( math::Sqrt( double( NUMBER_RESIDUES) / 100.0)));
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator that calculates quality measure of protein model to to coordinates stored in the aa sequences
    //! @param PROTEIN_MODEL the model to be considered
    //! @return the result of the quality measure this object was constructed with
    double Quality::operator()( const ProteinModel &PROTEIN_MODEL) const
    {
      return Calculate( m_Measure, PROTEIN_MODEL, m_AtomTypes);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read ProteinModel from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &Quality::Read( std::istream &ISTREAM)
    {
      // read member
      io::Serialize::Read( m_Measure, ISTREAM);
      io::Serialize::Read( m_AtomTypes, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write ProteinModel to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT nr of indentations
    //! @return output stream which was written to
    std::ostream &Quality::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write member
      io::Serialize::Write( m_Measure, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_AtomTypes, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief returns the pair of current and original coordinates for specified atomtypes for each amino acid that are in SSEs for the given chain
    //! @param CHAIN chain of interest
    //! @param NATIVE the native chain
    //! @param ATOM_TYPES Set of atom types to be returned
    //! @return the pair of current and original coordinates for specified atomtypes for each amino acid that are in SSEs for the given chain
    storage::Pair< util::SiPtrVector< const linal::Vector3D>, util::SiPtrVector< const linal::Vector3D> >
    Quality::GetSSEAndSequenceCoordinates
    (
      const Chain &CHAIN,
      const Chain &NATIVE,
      const storage::Set< biol::AtomType> &ATOM_TYPES
    )
    {
      //instantiate Pair of util::SiPtrVector of coordinates to model and according sequence coordinates
      storage::Pair< util::SiPtrVector< const linal::Vector3D>, util::SiPtrVector< const linal::Vector3D> >
        model_sequence_coordinates;

      //instantiate an iterator to sequence of amino acids
      biol::AASequence::const_iterator sequence_itr( NATIVE.GetSequence()->Begin());
      const biol::AASequence::const_iterator sequence_itr_end( NATIVE.GetSequence()->End());

      // instantiate secondary structure element iterator to the elements of the protein model
      for
      (
        storage::Set< util::ShPtr< SSE>, SSELessThanNoOverlap>::const_iterator
          sse_itr( CHAIN.GetData().Begin()), sse_itr_end( CHAIN.GetData().End());
        sse_itr != sse_itr_end && sequence_itr != sequence_itr_end;
        ++sse_itr
      )
      {
        //iterator to the first amino acids of current sselement
        biol::AASequence::const_iterator aa_itr( ( *sse_itr)->Begin());
        const biol::AASequence::const_iterator aa_itr_end( ( *sse_itr)->End());

        // increase sequence_itr till it matches the first aminoacid of the current SSE
        while
        (
          sequence_itr != sequence_itr_end &&
          aa_itr != aa_itr_end &&
          (
            ( *sequence_itr)->GetSeqID() < ( *aa_itr)->GetSeqID() ||
            (
              ( *sequence_itr)->GetSeqID() == ( *aa_itr)->GetSeqID() &&
              ( *sequence_itr)->GetPdbICode() < ( *aa_itr)->GetPdbICode()
            )
          )
        )
        {
          sequence_itr++;
        }

        // collect coordinates of common amino acids
        while
        (
          sequence_itr != sequence_itr_end &&
          aa_itr != aa_itr_end &&
          ( *sequence_itr)->GetSeqID() == ( *aa_itr)->GetSeqID() &&
          ( *sequence_itr)->GetPdbICode() == ( *aa_itr)->GetPdbICode()
        )
        {
          // initialize vectors to store coordinates for given atom types for original and current coordinates
          const util::SiPtrVector< const linal::Vector3D>
            model_coordinates( ( *aa_itr)->GetAtomCoordinates( ATOM_TYPES));
          const util::SiPtrVector< const linal::Vector3D> seq_coordinates
            ( ( *sequence_itr)->GetAtomCoordinates( ATOM_TYPES));

          //only if coordinates of current AA is defined
          if( coord::AreDefinedCoordinates( model_coordinates) && coord::AreDefinedCoordinates( seq_coordinates))
          {
            //store model coordinates
            model_sequence_coordinates.First().Append( model_coordinates);
            model_sequence_coordinates.Second().Append( seq_coordinates);
          }

          //goto next amino acid in sselement and sequence
          aa_itr++;
          sequence_itr++;
        }
      }

      // end
      return model_sequence_coordinates;
    }

    //! @brief returns the pair of current and original coordinates for specified atomtypes for each amino acid that are in SSEs for all the given chains
    //! @param MODEL Protein model of interest
    //! @param ATOM_TYPES Set of atom types to be returned
    //! @return the pair of current and original coordinates for specified atomtypes for each amino acid that are in SSEs for all the given chains
    storage::Pair< util::SiPtrVector< const linal::Vector3D>, util::SiPtrVector< const linal::Vector3D> >
    Quality::GetChainAndSequenceCoordinates
    (
      const ProteinModel &MODEL,
      const storage::Set< biol::AtomType> &ATOM_TYPES
    )
    {
      //instantiate Pair of util::SiPtrVector of coordinates to model and according sequence coordinates
      storage::Pair< util::SiPtrVector< const linal::Vector3D>, util::SiPtrVector< const linal::Vector3D> >
        model_sequence_coordinates;
      const util::SiPtr< const ProteinModel> template_model
      (
        MODEL.GetProteinModelData()->GetData( ProteinModelData::e_NativeModel)
      );

      // check that atom types are given
      if( ATOM_TYPES.IsEmpty())
      {
        BCL_MessageCrt( "no atom types given for rmsd calcualtion");
        return model_sequence_coordinates;
      }

      if( !template_model.IsDefined())
      {
        BCL_MessageCrt( "No native model given!");
        return model_sequence_coordinates;
      }

      // iterate over each chain
      for
      (
        util::ShPtrVector< Chain>::const_iterator
          chain_itr( MODEL.GetChains().Begin()), chain_itr_end( MODEL.GetChains().End());
        chain_itr != chain_itr_end; ++chain_itr
      )
      {
        util::ShPtr< Chain> sp_templ_chain( template_model->GetChain( ( *chain_itr)->GetChainID()));
        if( !sp_templ_chain.IsDefined())
        {
          BCL_MessageCrt( "Template model lacked chain " + std::string( size_t( 1), ( *chain_itr)->GetChainID()));
          return model_sequence_coordinates;
        }

        // get the current and original coordinates for the common amino acids
        storage::Pair< util::SiPtrVector< const linal::Vector3D>, util::SiPtrVector< const linal::Vector3D> >
          tmp_model_sequence_coordinates( GetSSEAndSequenceCoordinates( **chain_itr, *sp_templ_chain, ATOM_TYPES));

        // append to already stored coordinates
        model_sequence_coordinates.First().Append( tmp_model_sequence_coordinates.First());
        model_sequence_coordinates.Second().Append( tmp_model_sequence_coordinates.Second());
      }

      // end
      return model_sequence_coordinates;
    }

    //! @brief GetAllDesiredAtomCoordinates puts all the desired atoms for both proteins into a StorageVectorND
    //! @param COMMON_AAS VectorND which holds the amino acids which are common to both proteins
    //! @param DESIRED_ATOMS Set which holds the atoms which are desired to be included in the calculation
    //! @return VectorND< 2, util::SiPtrVector< const linal::Vector3D> > which has desired atom coordinates
    storage::VectorND< 2, util::SiPtrVector< const linal::Vector3D> > Quality::GetAllDesiredAtomCoordinates
    (
      const storage::VectorND< 2, util::SiPtrList< const biol::AABase> > &COMMON_AAS,
      const storage::Set< biol::AtomType> &DESIRED_ATOMS
    )
    {
      // initialize atom coordinates to return
      storage::VectorND< 2, util::SiPtrVector< const linal::Vector3D> > all_atom_coordinates;

      // get the coordinates of the atoms that are of interest
      for
      (
        util::SiPtrList< const biol::AABase>::const_iterator
          itr_a( COMMON_AAS.First().Begin()), itr_a_end( COMMON_AAS.First().End()),
          itr_b( COMMON_AAS.Second().Begin()), itr_b_end( COMMON_AAS.Second().End());
        itr_a != itr_a_end && itr_b != itr_b_end;
        ++itr_a, ++itr_b
      )
      {
        // add coordinates to the end of appropriate SiPtrVectors of "COMMON_AAS"
        all_atom_coordinates.First().Append( ( *itr_a)->GetAtomCoordinates( DESIRED_ATOMS));
        all_atom_coordinates.Second().Append( ( *itr_b)->GetAtomCoordinates( DESIRED_ATOMS));
      }

      // end
      return all_atom_coordinates;
    }

    //! @brief creates alignments from a single ProteinModel using the SSE residues
    //! @param PROTEIN_MODEL ProteinModel of interest
    //! @return alignments created from the ProteinModel using the SSE residues
    storage::Map< char, util::ShPtr< align::AlignmentInterface< biol::AABase> > >
    Quality::CreateAlignmentFromProteinModelSSEs( const ProteinModel &PROTEIN_MODEL)
    {
      // initialize map
      storage::Map< char, util::ShPtr< align::AlignmentInterface< biol::AABase> > > alignments_from_sses;

      // iterate over chains
      for
      (
        util::ShPtrVector< Chain>::const_iterator
          chain_itr( PROTEIN_MODEL.GetChains().Begin()), chain_itr_end( PROTEIN_MODEL.GetChains().End());
        chain_itr != chain_itr_end; ++chain_itr
      )
      {
        // initialize AA vector
        util::ShPtrVector< biol::AABase> residues;

        // iterate over the SSEs
        for
        (
          storage::Set< util::ShPtr< SSE>, SSELessThanNoOverlap>::const_iterator
            sse_itr( ( *chain_itr)->GetData().Begin()), sse_itr_end( ( *chain_itr)->GetData().End());
          sse_itr != sse_itr_end; ++sse_itr
        )
        {
          // add the residues to the vector
          residues.Append( ( *sse_itr)->GetData());
        }

        // create a new sequence
        const util::ShPtr< align::Sequence< biol::AABase> > sp_sequence( new align::Sequence< biol::AABase>( residues, ( *chain_itr)->GetSequence()->GetSequenceId()));

        // add to the map
        alignments_from_sses[ ( *chain_itr)->GetChainID()] =
          util::ShPtr< align::AlignmentInterface< biol::AABase> >( new align::AlignmentLeaf< biol::AABase>( sp_sequence));
      }

      // create and return the alignment
      return alignments_from_sses;
    }

    //! @brief creates alignments from a single ProteinModel using specific residues identified by collector
    //! @param PROTEIN_MODEL ProteinModel of interest
    //! @param COLLECTOR the collector that has the specific residues that will be used
    //! @return alignments created from the ProteinModel using the specified residues
    storage::Map< char, util::ShPtr< align::AlignmentInterface< biol::AABase> > >
    Quality::CreateAlignmentFromCollectorAASpecified
    (
      const ProteinModel &PROTEIN_MODEL, const CollectorAASpecified &COLLECTOR
    )
    {
      // get the residues of interest
      const util::SiPtrList< const biol::AABase> collected_aas( COLLECTOR.Collect( PROTEIN_MODEL));

      // will hold the residues for each chain
      storage::Map< char, util::ShPtrVector< biol::AABase> > chain_sequence;

      // iterate through the residues to build up the list per chain
      for
      (
        util::SiPtrList< const biol::AABase>::const_iterator
          aa_itr( collected_aas.Begin()), aa_itr_end( collected_aas.End());
        aa_itr != aa_itr_end;
        ++aa_itr
      )
      {
        // get reference to current aa list corresponding to the chain of the current residue
        util::ShPtrVector< biol::AABase> &residue_vector( chain_sequence[ ( *aa_itr)->GetChainID()]);

        // make shptr to cloned residue
        const util::ShPtr< biol::AABase> new_aa_ptr( ( *aa_itr)->Clone());

        residue_vector.PushBack( new_aa_ptr);

        if( residue_vector.GetSize() > 1)
        {
          BCL_Assert
          (
            new_aa_ptr->GetSeqID() > residue_vector( residue_vector.GetSize() - 2)->GetSeqID(),
            "the residues are being collected in non numerically increasing order. please make sure your desired "
            "residues to collect are listed in increasing numerical order."
          );
        }
      }

      // initialize map to hold chain and alignment
      storage::Map< char, util::ShPtr< align::AlignmentInterface< biol::AABase> > > alignments_from_sses;

      // iterate through the chain and sequence map
      for
      (
        storage::Map< char, util::ShPtrVector< biol::AABase> >::const_iterator
          chain_itr( chain_sequence.Begin()), chain_itr_end( chain_sequence.End());
        chain_itr != chain_itr_end;
        ++chain_itr
      )
      {
        // create a new sequence
        const util::ShPtr< align::Sequence< biol::AABase> > sp_sequence
        (
          new align::Sequence< biol::AABase>( chain_itr->second, std::string( &chain_itr->first, 1))
        );

        // add to the map
        alignments_from_sses[ chain_itr->first] =
          util::ShPtr< align::AlignmentInterface< biol::AABase> >( new align::AlignmentLeaf< biol::AABase>( sp_sequence));
      }

      return alignments_from_sses;
    }

    //! @brief creates alignments from a single ProteinModel using the chain sequence residues
    //! @param PROTEIN_MODEL ProteinModel of interest
    //! @return alignments created from the ProteinModel using the chain sequence residues
    storage::Map< char, util::ShPtr< align::AlignmentInterface< biol::AABase> > >
    Quality::CreateAlignmentFromProteinModelSequences( const ProteinModel &PROTEIN_MODEL)
    {
      // initialize map
      storage::Map< char, util::ShPtr< align::AlignmentInterface< biol::AABase> > > alignments_from_sequences;

      // iterate over chains
      for
      (
        util::ShPtrVector< Chain>::const_iterator
          chain_itr( PROTEIN_MODEL.GetChains().Begin()), chain_itr_end( PROTEIN_MODEL.GetChains().End());
        chain_itr != chain_itr_end; ++chain_itr
      )
      {
        // create an alignment from the chain sequence
        const util::ShPtr< biol::AASequence> &sp_sequence( ( *chain_itr)->GetSequence());
        util::ShPtr< align::AlignmentLeaf< biol::AABase> > sp_leaf
        (
          new align::AlignmentLeaf< biol::AABase>( sp_sequence)
        );

        // add the alignment to the map
        alignments_from_sequences[ ( *chain_itr)->GetChainID()] = sp_leaf;
      }

      // end
      return alignments_from_sequences;
    }

    //! @brief creates alignments of chain sequences to SSE residues for a given ProteinModel
    //! @param PROTEIN_MODEL ProteinModel of interest
    //! @return alignments of chain sequences to SSE residues for a given ProteinModel
    storage::Map< char, util::ShPtr< align::AlignmentInterface< biol::AABase> > >
    Quality::CreateAlignmentProteinModelSequencesToSSEs( const ProteinModel &PROTEIN_MODEL)
    {
      // get the alignments of sequences
      storage::Map< char, util::ShPtr< align::AlignmentInterface< biol::AABase> > > align_seq
      (
        Quality::CreateAlignmentFromProteinModelSequences( PROTEIN_MODEL)
      );

      // get the alignments of SSEs
      storage::Map< char, util::ShPtr< align::AlignmentInterface< biol::AABase> > > align_sses
      (
        Quality::CreateAlignmentFromProteinModelSSEs( PROTEIN_MODEL)
      );

      // return the combined alignment
      return AlignSingleAlignments( align_seq, align_sses);
    }

    //! @brief creates alignments of sequences in SSEs of two protein models
    //! @param PROTEIN_MODEL_A first ProteinModel of interest
    //! @param PROTEIN_MODEL_B second ProteinModel of interest
    //! @return alignments of sequences in SSEs of two protein models
    storage::Map< char, util::ShPtr< align::AlignmentInterface< biol::AABase> > >
    Quality::CreateAlignmentProteinModels
    (
      const ProteinModel &PROTEIN_MODEL_A,
      const ProteinModel &PROTEIN_MODEL_B
    )
    {
      // get the alignment of model a
      storage::Map< char, util::ShPtr< align::AlignmentInterface< biol::AABase> > > align_a
      (
        Quality::CreateAlignmentFromProteinModelSSEs( PROTEIN_MODEL_A)
      );

      // get the alignments of model b
      storage::Map< char, util::ShPtr< align::AlignmentInterface< biol::AABase> > > align_b
      (
        Quality::CreateAlignmentFromProteinModelSSEs( PROTEIN_MODEL_B)
      );

      // return the combined alignment
      return AlignSingleAlignments( align_a, align_b, biol::AACompareData());
    }

    //! @brief aligns two single alignments (each containing 1 sequence) using the given comparison
    //! @param ALIGNMENT_A first alignment
    //! @param ALIGNMENT_B second alignment
    //! @param COMPARISON comparison object for two AABase
    //! @return combined alignments
    storage::Map< char, util::ShPtr< align::AlignmentInterface< biol::AABase> > >
    Quality::AlignSingleAlignments
    (
      storage::Map< char, util::ShPtr< align::AlignmentInterface< biol::AABase> > > &ALIGNMENT_A,
      storage::Map< char, util::ShPtr< align::AlignmentInterface< biol::AABase> > > &ALIGNMENT_B,
      const util::BinaryFunctionInterface< biol::AABase, biol::AABase, bool> &COMPARISON
    )
    {
      function::BinarySum< const biol::AABase, const biol::AABase, double> assign_score;
      assign_score += **score::GetAAAssignments().e_IDENTITY;

      // Based on Elizabeth Dong's Weight tables -0.672 -0.907 -0.635  0.00
      // fold recognition scores ENCLOSED_SINGLE_GAP, ENCLOSED_MULTIPLE_GAP, BOUNDARY_SINGLE_GAP, BOUNDARY_MULTIPLE_GAP
      score::AssignmentWithGap< biol::AABase> assign_gap_score
      (
        util::CloneToShPtr( assign_score),
        -0.672, // ENCLOSED_SINGLE_GAP
        -0.907, // ENCLOSED_MULTIPLE_GAP
        -0.635, // BOUNDARY_SINGLE_GAP
        0.0 // BOUNDARY_MULTIPLE_GAP
      );
      // just identity score
      score::AssignmentWithGap< biol::AABase> identity_score
      (
        util::CloneToShPtr( assign_score),
        0.0, // ENCLOSED_SINGLE_GAP
        0.0, // ENCLOSED_MULTIPLE_GAP
        0.0, // BOUNDARY_SINGLE_GAP
        0.0 // BOUNDARY_MULTIPLE_GAP
      );

      // initialize map
      storage::Map< char, util::ShPtr< align::AlignmentInterface< biol::AABase> > > sequences;

      // iterate over the model a map
      for
      (
        storage::Map< char, util::ShPtr< align::AlignmentInterface< biol::AABase> > >::iterator
          seq_itr( ALIGNMENT_A.Begin()), seq_itr_end( ALIGNMENT_A.End());
        seq_itr != seq_itr_end; ++seq_itr
      )
      {
        // try to find the corresponding sequence from the sse map
        storage::Map< char, util::ShPtr< align::AlignmentInterface< biol::AABase> > >::iterator find_itr
        (
          ALIGNMENT_B.Find( seq_itr->first)
        );

        // if the sequence was found
        if( find_itr != ALIGNMENT_B.End())
        {
          // add the alignment to the map
          sequences[ seq_itr->first] =
            util::ShPtr< align::AlignmentInterface< biol::AABase> >
            (
              biol::AlignByPdbID().AlignPair( seq_itr->second, find_itr->second).First().Clone()
            );
          if( sequences[ seq_itr->first]->Score( identity_score) < std::min( find_itr->second->GetSize(), seq_itr->second->GetSize()))
          {
            BCL_MessageVrb( "Recomputing alignments because PDB IDs of the native and the model did not match!")
            align::AlignerDP< biol::AABase> aligner;
            aligner.SetScoringFunction( assign_gap_score);
            sequences[ seq_itr->first] =
                        util::ShPtr< align::AlignmentInterface< biol::AABase> >
                        (
                          aligner.AlignPair( seq_itr->second, find_itr->second).First().Clone()
                        );
          }
        }
      }

      // end
      return sequences;
    }

    //! @brief prunes the given alignment by removing all assignments with gaps
    //! @param ALIGNMENT Alignment of two protein sequences
    //! @return pruned alignment where the gaps were removed
    align::AlignmentNode< biol::AABase> Quality::RemoveAssignmentsWithGapsFromAlignment
    (
      const align::AlignmentInterface< biol::AABase> &ALIGNMENT
    )
    {
      // initialize alignment
      align::AlignmentNode< biol::AABase> pruned_alignment( ALIGNMENT.GetChildAlignments());

      // iterate over the given alignment
      for
      (
        align::AlignmentLeaf< biol::AABase>::const_iterator
          itr( ALIGNMENT.GetAssignments().Begin()), itr_end( ALIGNMENT.GetAssignments().End());
        itr != itr_end; ++itr
      )
      {
        // if all the members of the assignment are defined
        if( ( *itr)->GetMembers().IsDefined())
        {
          pruned_alignment.Append( *itr);
        }
      }

      // end
      return pruned_alignment;
    }

    //! @brief prunes the given alignment making sure each template amino acid has defined coordinates for given atom types
    //! @param ALIGNMENT Alignment of two protein sequences
    //! @param ATOM_TYPES Set of atom types for which the coordinates need to be defined
    //! @return pruned alignment where the residues where the template sequence has undefined coordinates were removed
    align::AlignmentNode< biol::AABase> Quality::RemoveUndefinedTemplateAminoAcidsFromAlignment
    (
      const align::AlignmentInterface< biol::AABase> &ALIGNMENT,
      const storage::Set< biol::AtomType> &ATOM_TYPES
    )
    {
      // initialize alignment
      align::AlignmentNode< biol::AABase> pruned_alignment( ALIGNMENT.GetChildAlignments());

      // iterate over the given alignment
      for
      (
        align::AlignmentLeaf< biol::AABase>::const_iterator
          itr( ALIGNMENT.GetAssignments().Begin()), itr_end( ALIGNMENT.GetAssignments().End());
        itr != itr_end; ++itr
      )
      {
        // create a pointer on the template (first sequence) amino acid for this position
        const util::SiPtr< const biol::AABase> aa_ptr( ( *itr)->GetMembers().FirstElement());

        // if the amino acid is not defined, skip to next assignment
        if( !aa_ptr.IsDefined())
        {
          continue;
        }

        // get the coordinates for this amino acid
        const util::SiPtrVector< const linal::Vector3D> coords( aa_ptr->GetAtomCoordinates( ATOM_TYPES));

        // if coordinates are defined
        if( !coords.IsEmpty() && coord::AreDefinedCoordinates( coords))
        {
          // add this assignment to the pruned alignment created
          pruned_alignment.Append( *itr);
        }
      }

      if( ALIGNMENT.GetSize() != pruned_alignment.GetSize())
      {
        BCL_MessageStd
        (
          "Pruned alignment from " + util::Format()( ALIGNMENT.GetSize()) + " to "
          + util::Format()( pruned_alignment.GetSize()) + " assignments by removing undefined coordinates."
        );
      }
      // end
      return pruned_alignment;
    }

    //! @brief prunes the given alignment making sure each amino acid has defined coordinates for given atom types
    //! @param ALIGNMENT Alignment of two protein sequences
    //! @param ATOM_TYPES Set of atom types for which the coordinates need to be defined
    //! @return pruned alignment where the residues with undefined coordinates were removed
    align::AlignmentNode< biol::AABase> Quality::RemoveUndefinedAminoAcidsFromAlignment
    (
      const align::AlignmentInterface< biol::AABase> &ALIGNMENT,
      const storage::Set< biol::AtomType> &ATOM_TYPES
    )
    {
      // initialize alignment
      align::AlignmentNode< biol::AABase> pruned_alignment( ALIGNMENT.GetChildAlignments());

      // iterate over the given alignment
      for
      (
        align::AlignmentLeaf< biol::AABase>::const_iterator
          itr( ALIGNMENT.GetAssignments().Begin()), itr_end( ALIGNMENT.GetAssignments().End());
        itr != itr_end; ++itr
      )
      {
        // if all the members of the assignment are defined
        bool coordinates_defined( true);

        // iterate over all amino acids to check if the atom coordinates of interest are defined
        for
        (
          util::SiPtrList< const biol::AABase>::const_iterator
            itr_mem( ( *itr)->GetMembers().Begin()), itr_mem_end( ( *itr)->GetMembers().End());
          itr_mem != itr_mem_end && coordinates_defined;
          ++itr_mem
        )
        {
          // if amino acid is not defined
          if( !itr_mem->IsDefined())
          {
            // set corodinates_defined to false and break
            coordinates_defined = false;
            break;
          }

          // get the coordinates
          const util::SiPtrVector< const linal::Vector3D> coords( ( *itr_mem)->GetAtomCoordinates( ATOM_TYPES));

          coordinates_defined &= ( !coords.IsEmpty() && coord::AreDefinedCoordinates( coords));
        }

        // insert if coordinates are defined
        if( coordinates_defined)
        {
          pruned_alignment.Append( *itr);
        }
      }

      // end
      return pruned_alignment;
    }

    //! @brief given an amino acid alignment and atom types generate two siptr vector of corresponding coordinates
    //!        for quality calculations, and the number of residues involved (e.g. skipping undefined residues
    //! @param ALIGNMENT Alignment of two protein sequences
    //! @param ATOM_TYPES Set of atom types of interest
    //! @return two siptr vector of corresponding coordinates, and the number of residues in the alignment
    storage::Pair< storage::VectorND< 2, util::SiPtrVector< const linal::Vector3D> >, size_t>
    Quality::CoordinatesAndResidueCountFromAlignment
    (
      const align::AlignmentInterface< biol::AABase> &ALIGNMENT,
      const storage::Set< biol::AtomType> &ATOM_TYPES
    )
    {
      storage::Pair< storage::VectorND< 2, util::SiPtrVector< const linal::Vector3D> >, size_t>
      coordinates_pair_and_nres
      (
        storage::VectorND< 2, util::SiPtrVector< const linal::Vector3D> >(),
        size_t( 0)
      );
      size_t n_residues( 0);
      storage::VectorND< 2, util::SiPtrVector< const linal::Vector3D> > &coordinates_pair
      (
        coordinates_pair_and_nres.First()
      );

      if( ALIGNMENT.GetDepth() < 2)
      {
        return coordinates_pair_and_nres;
      }

      // allocate enough space
      coordinates_pair.First().AllocateMemory( ATOM_TYPES.GetSize() * ALIGNMENT.GetSize());
      coordinates_pair.Second().AllocateMemory( ATOM_TYPES.GetSize() * ALIGNMENT.GetSize());

      // iterate over the given alignment
      for
      (
        align::AlignmentLeaf< biol::AABase>::const_iterator
          itr( ALIGNMENT.GetAssignments().Begin()), itr_end( ALIGNMENT.GetAssignments().End());
        itr != itr_end; ++itr
      )
      {
        const util::SiPtrList< const biol::AABase> &members( ( *itr)->GetMembers());
        const util::SiPtr< const biol::AABase> &member_a( *members.Begin()); // first member
        const util::SiPtr< const biol::AABase> &member_b( *++members.Begin()); // second member

        // check that amino acids is not aligned against a gap
        if( !member_a.IsDefined() || !member_b.IsDefined())
        {
          continue;
        }

        // get the coordinates for the given atom types from the amino acids
        const util::SiPtrVector< const linal::Vector3D> coords_a( member_a->GetAtomCoordinates( ATOM_TYPES));
        const util::SiPtrVector< const linal::Vector3D> coords_b( member_b->GetAtomCoordinates( ATOM_TYPES));

        // check that the number of coordinates agrees and both are defined
        if
        (
             coords_a.GetSize() != coords_b.GetSize()
          || !coord::AreDefinedCoordinates( coords_a)
          || !coord::AreDefinedCoordinates( coords_b)
        )
        {
          continue;
        }

        // insert the coordinates
        ++n_residues;
        coordinates_pair.First().Append( coords_a);
        coordinates_pair.Second().Append( coords_b);
      }

      coordinates_pair_and_nres.Second() = n_residues;

      // end
      return coordinates_pair_and_nres;
    }

    //! @brief given an amino acid alignment and atom types generate two siptr vector of corresponding coordinates
    //!        for quality calculations
    //! @param ALIGNMENT Alignment of two protein sequences
    //! @param ATOM_TYPES Set of atom types of interest
    //! @return two siptr vector of corresponding coordinates
    storage::VectorND< 2, util::SiPtrVector< const linal::Vector3D> > Quality::CoordinatesFromAlignment
    (
      const align::AlignmentInterface< biol::AABase> &ALIGNMENT,
      const storage::Set< biol::AtomType> &ATOM_TYPES
    )
    {
      return CoordinatesAndResidueCountFromAlignment( ALIGNMENT, ATOM_TYPES).First();
    }

    //! @brief given an amino acid alignments and atom types generate two siptr vector of corresponding coordinates
    //!        for quality calculations
    //! @param ALIGNMENTS Alignments of two protein sequences
    //! @param ATOM_TYPES Set of atom types of interest
    //! @return two siptr vector of corresponding coordinates
    storage::VectorND< 2, util::SiPtrVector< const linal::Vector3D> > Quality::CoordinatesFromAlignments
    (
      const storage::Map< char, util::ShPtr< align::AlignmentInterface< biol::AABase> > > &ALIGNMENTS,
      const storage::Set< biol::AtomType> &ATOM_TYPES
    )
    {
      // initialize coordinates vector
      storage::VectorND< 2, util::SiPtrVector< const linal::Vector3D> > coords;

      // iterate through the alignment map
      for
      (
        storage::Map< char, util::ShPtr< align::AlignmentInterface< biol::AABase> > >::const_iterator
          align_itr( ALIGNMENTS.Begin()), align_itr_end( ALIGNMENTS.End());
        align_itr != align_itr_end; ++align_itr
      )
      {
        // get the common coordinates
        const storage::VectorND< 2, util::SiPtrVector< const linal::Vector3D> > this_coordinates_pair
        (
          CoordinatesFromAlignment( *( align_itr->second), ATOM_TYPES)
        );

        // update the total coordinates
        coords.First().Append( this_coordinates_pair.First());
        coords.Second().Append( this_coordinates_pair.Second());
      }

      // end
      return coords;
    }

    //! @brief uses a vector of filenames to read in alignments for chains - chain is last character before extension
    //! @param HANDLER the type of handler to read in the alignments
    //! @param ALIGN_FILES the files containing the alignments for each chain
    //! @return map with a character for the chain and a corresponding alignment for that chain
    storage::Map< char, util::ShPtr< align::AlignmentInterface< biol::AABase> > > Quality::GetChainAlignments
    (
      const align::HandlerInterface< biol::AABase> &HANDLER, const storage::Vector< io::DirectoryEntry> &ALIGN_FILES
    )
    {
      // to hold for each chain an alignment
      storage::Map< char, util::ShPtr< align::AlignmentInterface< biol::AABase> > > chain_align;

      // iterate over the alignment files
      for
      (
        storage::Vector< io::DirectoryEntry>::const_iterator align_itr( ALIGN_FILES.Begin()),
          align_itr_end( ALIGN_FILES.End());
        align_itr != align_itr_end; ++align_itr
      )
      {
        // get the filename without extensions
        const std::string extensionless_file( io::File::RemoveFullExtension( align_itr->GetName()));

        // get the chain id - last character of extension-less filename
        const char chain_id( extensionless_file[ extensionless_file.length() - 1]);

        // read in the alignment
        io::IFStream read;
        io::File::MustOpenIFStream( read, align_itr->GetFullName());
        chain_align[ chain_id] = HANDLER.ReadAlignment( read, biol::AASequence(), chain_id);
        io::File::CloseClearFStream( read);
      }

      return chain_align;
    }

  } // namespace assemble
} // namespace bcl

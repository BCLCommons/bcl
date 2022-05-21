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
#include "fold/bcl_fold_add_parabolic_loops.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_atom.h"
#include "opti/bcl_opti_approximator_root_regula_falsi.h"
#include "opti/bcl_opti_criterion_combine.h"
#include "opti/bcl_opti_criterion_convergence_argument.h"
#include "opti/bcl_opti_criterion_number_iterations.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class NormOptimization
    //! @brief takes a normalization factor as input and calculates the parabolic arc length - desired arc length
    //! @detail the arc length is computed using the formula:
    //!     L = 1/4a[ ln((d2+s2)/(d1+s1)) + d2s2 -d1s1],
    //!     l1 = 0.0              l2 = 1.0
    //!      a = norm_factor       b = -(norm_factor + SSE Distance)
    //!     d1 = 2al1 + b         d2 = 2al2 + b
    //!     s1 = Sqrt( 1 + d1^2)  s2 = Sqrt( 1 + d2^2)
    //!AddParabolicLoops::
    //! @see @link example_fold_add_parabolic_loops.cpp @endlink
    //! @author putnamdk, mendenjl
    //! @date Aug 28, 2012
    //!
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    class BCL_API NormOptimization :
      public math::FunctionInterfaceSerializable< double, double>
    {

    private:

    //////////
    // data //
    //////////

      double m_DesiredDistance; //!< The desired parabolic length
      double m_SSEDistance;     //!< The linear distance between sse's

    public:

    //////////
    // constructors //
    //////////

      //! @brief Clone function
      //! @return pointer to new NormOptimization
      NormOptimization *Clone() const
      {
        return new NormOptimization( *this);
      }

      //! @brief constructor from a linear distance between sse's
      //! @param DESIRED_DISTANCE double which is the desired length of the parabolic path between two SSE's
      //! @param SSE_DISTANCE double which is the linear distance between two SSE's
      NormOptimization( const double DESIRED_DISTANCE, const double SSE_DISTANCE) :
        m_DesiredDistance( DESIRED_DISTANCE),
        m_SSEDistance( SSE_DISTANCE)
      {
      }

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const
      {
        // Get BCL Standardized Class name
        return GetStaticClassName( *this);
      }

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM)
      {
        // return the stream
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        // return the stream
        return OSTREAM;
      }

      //! @brief takes a normalization factor as input and calculates the parabolic arc length - desired arc length
      //! @detail the arc length is computed using the formula:
      //!     L = 1/4a[ ln((d2+s2)/(d1+s1)) + d2s2 -d1s1],
      //!     l1 = 0.0              l2 = 1.0
      //!      a = norm_factor       b = -(norm_factor + SSE Distance)
      //!     d1 = 2al1 + b         d2 = 2al2 + b
      //!     s1 = Sqrt( 1 + d1^2)  s2 = Sqrt( 1 + d2^2)
      //! @param ARGUMENT a given normalization factor
      //! @return the parabolic arc length - the desired arc length for a given normalization factor
      double operator()( const double &ARGUMENT) const
      {
        const double a( ARGUMENT);
        const double b( -( ARGUMENT + m_SSEDistance));

        const double partial_derivative( 2 * a + b);

        const double intermediate_1( math::Sqrt( 1 + math::Sqr( b)));
        const double intermediate_2( math::Sqrt( 1 + math::Sqr( partial_derivative)));

        const double constant( 1.0 / ( 4.0 * a));
        const double log_part( std::log( ( partial_derivative + intermediate_2) / ( b + intermediate_1)));
        const double parabolic_length( constant * ( log_part + ( partial_derivative * intermediate_2 - b * intermediate_1)));

        return parabolic_length - m_DesiredDistance;
      }
    };

  //////////
  // data //
  //////////

    //! single instance of that class
    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> AddParabolicLoops::s_Instance
    (
      GetObjectInstances().AddInstance( new AddParabolicLoops( true))
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief Constructor that takes a bool
    //! @param LOOPS bool value to represent loops that are not present in the protein model
    AddParabolicLoops::AddParabolicLoops( const bool &ANALYTIC_NORM_FACTOR, const bool &APPROXIMATE_CB) :
      m_DetermineAnalyticNormFactor( ANALYTIC_NORM_FACTOR),
      m_ApproximateCB( APPROXIMATE_CB)
    {
    }

    //! @brief Clone function
    //! @return pointer add parabolic loops
    AddParabolicLoops *AddParabolicLoops::Clone() const
    {
      return new AddParabolicLoops( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AddParabolicLoops::GetClassIdentifier() const
    {
      // Get BCL Standardized Class name
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

  ///////////////
  // operators //
  ///////////////

    //! @brief TODO
    //! @param PROTEIN_MODEL - protein model
    //! @param FIRST_RESIDUE, LAST_RESIDUE, first and last residues for the loop
    //! @return returns the coordinates for each atom in the loop
    storage::Vector< storage::Pair< biol::AAType, linal::Vector3D> >
    AddParabolicLoops::GetLoopCoordinates
    (
      const assemble::ProteinModel &PROTEIN_MODEL
    ) const
    {
      storage::Vector< storage::Pair< biol::AAType, linal::Vector3D> > loop_coordinates_and_atom_types;

      // iterate over the chains
      for
      (
        util::ShPtrVector< assemble::Chain>::const_iterator chain_itr( PROTEIN_MODEL.GetChains().Begin()),
          chain_itr_end( PROTEIN_MODEL.GetChains().End());
        chain_itr != chain_itr_end; ++chain_itr
      )
      {
        // get all the structured sses
        const util::SiPtrVector< const assemble::SSE> structured_sses
        (
          ( *chain_itr)->GetSSEs( storage::Set< biol::SSType>( biol::GetSSTypes().HELIX, biol::GetSSTypes().STRAND))
        );

        // if structured_sses are empty
        if( structured_sses.IsEmpty())
        {
          // warn user and return empty vector
          BCL_MessageStd( "No structured SSEs found in protein model");
          return loop_coordinates_and_atom_types;
        }

        // get the first and last sse
        const util::SiPtr< const assemble::SSE> sp_first_sse( structured_sses.FirstElement());
        const util::SiPtr< const assemble::SSE> sp_last_sse( structured_sses.LastElement());

        int prev_seqid( sp_first_sse->GetLastAA()->GetSeqID());
        linal::Vector3D prev_sse_end( sp_first_sse->GetSSEGeometries().LastElement()->EndOfZ());
        linal::Vector3D prev_sse_ctr( sp_first_sse->GetSSEGeometries().LastElement()->GetCenter());

        // iterate over the SSEs
        for
        (
          storage::Set< util::ShPtr< assemble::SSE>, assemble::SSELessThanNoOverlap>::const_iterator
            sse_itr( ( *chain_itr)->GetData().Begin()), sse_itr_end( ( *chain_itr)->GetData().End());
          sse_itr != sse_itr_end; ++sse_itr
        )
        {
          // if estimating loops and this sse is not a coil
          if( ( *sse_itr)->GetType()->IsStructured())
          {
            // store the begin and end id
            const int begin_seqid( ( *sse_itr)->GetFirstAA()->GetSeqID());

            // if there is any region between this and the previous sse
            if( begin_seqid > prev_seqid + 1)
            {
              // get the sequence
              const util::ShPtr< biol::AASequence> &chain_sequence( ( *chain_itr)->GetSequence());
              biol::AASequence sub_sequence
              (
                chain_sequence->SubSequence( prev_seqid, begin_seqid - prev_seqid - 1)
              );

              // Calculate the length of the subsequence
              const size_t sub_sequence_length( sub_sequence.GetSize());

              //Calculate Vector
              const linal::Vector3D begin_coords( ( *sse_itr)->GetSSEGeometries().FirstElement()->BeginOfZ());

              // The coordinates will be approximated by parabolic loop building function
              // Calculate Vector pointing in the direction of the SSE by subtracting Initialize starting loop point
              const linal::Vector3D point_a( prev_sse_end);
              const linal::Vector3D center_a( prev_sse_ctr);
              const linal::Vector3D vector_a( point_a - center_a);

              // Get the direction of the vector pointing along the first SSE by diving all elements by the magnitude
              const linal::Vector3D dir_vector_a( vector_a / vector_a.Norm());

              // Initialize the end loop point
              // Calculate Vector pointing in the direction of the SSE by subtracting Initialize starting loop point
              const linal::Vector3D point_b( begin_coords);
              const linal::Vector3D center_b( ( *sse_itr)->GetSSEGeometries().FirstElement()->GetCenter());
              const linal::Vector3D vector_b( point_b - center_b);

              // Get the direction of the vector pointing along the second SSE by diving all elements by the magnitude
              const linal::Vector3D dir_vector_b( vector_b / vector_b.Norm());

              // Find the distance between points along the parabolic path.  This is a number between 0 and 1
              double distance = 1.0 / double( sub_sequence_length + 1);

              // Find the vector between Secondary Structure elements
              const linal::Vector3D vector_c( point_b - point_a);

              // Compute the distance between Secondary Structure elements
              const double sse_distance( vector_c.Norm());

              // Set the CB-CB spacing
              const double atom_distance( 3.2);

              // Compute the desired distance
              const double desired_distance( ( sub_sequence_length) * atom_distance);

              const double norm_factor
              (
                m_DetermineAnalyticNormFactor
                ? ComputeNormFactor( desired_distance, sse_distance)
                : ComputeNormalizationFactor( desired_distance, sse_distance, sub_sequence_length)
              );

              // Initialize the location of the point
              double location = 0.0;

              // Initialize the point to 0
              linal::Vector3D point( 0, 0, 0);

              BCL_MessageDbg( "Found loop with size: " + util::Format()( sub_sequence.GetData().GetSize()));

              // Iterate over the subsequence
              for
              (
                util::ShPtrVector< biol::AABase>::iterator sub_itr( sub_sequence.GetData().Begin()),
                  sub_itr_end( sub_sequence.GetData().End());
                sub_itr != sub_itr_end;
                ++sub_itr
              )
              {
                location += distance;
                point = ( 1.0 - location) * point_a + location * point_b +
                          norm_factor * location * ( 1.0 - location) *
                        ( ( 1.0 - location) * dir_vector_a + ( location * dir_vector_b));

                loop_coordinates_and_atom_types.PushBack
                (
                  storage::Pair< biol::AAType, linal::Vector3D>( ( *sub_itr)->GetType(), point)
                );
              }
            } // end loop region check

            // update previous variables
            prev_seqid = ( *sse_itr)->GetLastAA()->GetSeqID();
            prev_sse_end = ( *sse_itr)->GetSSEGeometries().LastElement()->EndOfZ();
            prev_sse_ctr = ( *sse_itr)->GetSSEGeometries().LastElement()->GetCenter();
          } // end loop estimation
        } // SSE iteration
      } // chain iteration

      return loop_coordinates_and_atom_types;
    }

    //! @brief takes a proteinModel as input and calculates an intensity using the debye formula
    //! @param PROTEIN_MODEL
    //! @return the intensity for a given q value
    assemble::ProteinModel AddParabolicLoops::operator()( const assemble::ProteinModel &PROTEIN_MODEL) const
    {
      util::ShPtr< assemble::ProteinModel> sp_tmp_protein( PROTEIN_MODEL.HardCopy());
      assemble::ProteinModel &tmp_protein( *sp_tmp_protein);

      storage::Vector< storage::Pair< biol::AAType, linal::Vector3D> > loop_coordinates_aa_types
      (
        GetLoopCoordinates( tmp_protein)
      );

      storage::Vector< storage::Pair< biol::AAType, linal::Vector3D> >::const_iterator
        itr_loop_coordinates_aa_types( loop_coordinates_aa_types.Begin()),
        itr_loop_coordinates_aa_types_end( loop_coordinates_aa_types.End());

      // iterate over the chains
      for
      (
        util::ShPtrVector< assemble::Chain>::iterator chain_itr( tmp_protein.GetChains().Begin()),
          chain_itr_end( tmp_protein.GetChains().End());
        chain_itr != chain_itr_end; ++chain_itr
      )
      {
        ( *chain_itr)->AddLoops( true, false);

        // get all the structured sses
        const util::SiPtrVector< const assemble::SSE> unstructured_sses
        (
          ( *chain_itr)->GetSSEs( biol::GetSSTypes().COIL)
        );

        // get the first and last residues' seq-ids
        const int first_seq_id( ( *chain_itr)->GetSequence()->GetData().FirstElement()->GetSeqID());
        const int last_seq_id( ( *chain_itr)->GetSequence()->GetData().LastElement()->GetSeqID());

        // iterate over the SSEs
        for
        (
          storage::Set< util::ShPtr< assemble::SSE>, assemble::SSELessThanNoOverlap>::const_iterator
            sse_itr( ( *chain_itr)->GetData().Begin()), sse_itr_end( ( *chain_itr)->GetData().End());
          sse_itr != sse_itr_end; ++sse_itr
        )
        {
          util::ShPtr< assemble::SSE> sse( *sse_itr);

          // if estimating loops and this sse is not a coil
          if( !sse->GetType()->IsStructured())
          {
            if( sse->GetLastAA()->GetSeqID() == last_seq_id || sse->GetFirstAA()->GetSeqID() == first_seq_id)
            {
              // skip terminal loops
              continue;
            }
            // Iterate over the subsequence
            BCL_MessageDbg( "Setting coordinates for loop with size: " + util::Format()( sse->GetData().GetSize()));
            math::RunningAverage< linal::Vector3D> loop_center;
            if( m_ApproximateCB)
            {
              size_t loop_size( sse->GetData().GetSize());
              auto itr_loop_coord_aa_types_cp( itr_loop_coordinates_aa_types);
              for( size_t pos( 0); pos < loop_size; ++pos, ++itr_loop_coord_aa_types_cp)
              {
                loop_center += itr_loop_coord_aa_types_cp->Second();
              }
            }

            for
            (
              util::ShPtrVector< biol::AABase>::iterator sub_itr( sse->GetData().Begin()),
                sub_itr_end( sse->GetData().End());
              sub_itr != sub_itr_end;
              ++sub_itr, ++itr_loop_coordinates_aa_types
            )
            {
              BCL_Assert
              (
                itr_loop_coordinates_aa_types != itr_loop_coordinates_aa_types_end,
                "Premature end to loop coordinates vector " + util::Format()( loop_coordinates_aa_types.GetSize())
                + " loop size " + util::Format()( sse->GetData().GetSize())
              );

              const linal::Vector3D &coordinate( itr_loop_coordinates_aa_types->Second());

              // create atoms at interpolated coordinates, set atoms in amino acid
              biol::Atom atom_n ( coordinate, biol::GetAtomTypes().N,  util::GetUndefinedSize_t(), -1.00);
              biol::Atom atom_ca( coordinate, biol::GetAtomTypes().CA, util::GetUndefinedSize_t(), -1.00);
              biol::Atom atom_c ( coordinate, biol::GetAtomTypes().C,  util::GetUndefinedSize_t(), -1.00);
              biol::Atom atom_o ( coordinate, biol::GetAtomTypes().O,  util::GetUndefinedSize_t(), -1.00);
              biol::Atom atom_cb( coordinate, biol::GetAtomTypes().CB, util::GetUndefinedSize_t(), -1.00);
              if( m_ApproximateCB)
              {
                atom_cb = biol::Atom( linal::UnitVector( loop_center, coordinate) * 1.4 + coordinate, biol::GetAtomTypes().CB, util::GetUndefinedSize_t(), -1.00);
              }
              ( **sub_itr).SetAtoms
              (
                util::SiPtrVector< const biol::Atom>::Create( atom_n, atom_ca, atom_c, atom_o, atom_cb)
              );
            }
          } // end loop estimation
        } // SSE iteration
      } // chain iteration

      BCL_Assert
      (
        itr_loop_coordinates_aa_types == itr_loop_coordinates_aa_types_end,
        "Should have used all the loop coordinates, but omitted "
        + util::Format()( std::distance( itr_loop_coordinates_aa_types, itr_loop_coordinates_aa_types_end))
        + " of them"
      );
      return tmp_protein;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &AddParabolicLoops::Read( std::istream &ISTREAM)
    {
      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &AddParabolicLoops::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return Computed Normalization factor for parabolic function - determines parabolic height.
    //! @param DESIRED_DISTANCE - this is how long you want the parabola to be
    //! @param SSE_DISTANCE - this is the distance between two sse's
    //! @return double normalization factor
    double AddParabolicLoops::ComputeNormFactor( const double &DESIRED_DISTANCE, const double &SSE_DISTANCE) const
    {
      // Set the distance from 0 the function will find
      const double root_tolerance( 0.0001);

      // Set Left point of interval, must be > 0
      const double border_left( root_tolerance);

      // Set Right point of interval, must be sufficient to contain x-intercept
      const double border_right( 5.0);

      // build the objective function
      util::ShPtr< math::FunctionInterfaceSerializable< double, double> >
        sp_function( new NormOptimization( DESIRED_DISTANCE, SSE_DISTANCE));

      // combine the termination criteria
      opti::CriterionCombine< double, double> criterion_combine;
      criterion_combine.InsertCriteria( opti::CriterionNumberIterations< double, double>( 200));
      criterion_combine.InsertCriteria( opti::CriterionConvergenceArgument< double, double>( 1, root_tolerance));

      // create regula falsi approximator
      opti::ApproximatorRootRegulaFalsi< double, double> approximator
      (
        *sp_function, criterion_combine, border_left, border_right
      );

      // approximate
      approximator.Approximate();

      // return the x-value that approximates the function value to 0 with the set tolerance
      return approximator.GetTracker().GetBest()->First();
    }

    double AddParabolicLoops::ComputeNormalizationFactor( const double &DESIRED_DISTANCE, const double &SSE_DISTANCE, const double &NAA) const
    {
      if( DESIRED_DISTANCE > SSE_DISTANCE)
      {
        const double height( 0.5 * math::Sqrt( math::Sqr( DESIRED_DISTANCE) - math::Sqr( SSE_DISTANCE)));
        return height;
      }
      return 0.0;
    }
  } // namespace fold
} // namespace bcl

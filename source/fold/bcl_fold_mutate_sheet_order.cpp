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
#include "fold/bcl_fold_mutate_sheet_order.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_domain.h"
#include "coord/bcl_coord_move_rotate_defined.h"
#include "io/bcl_io_serialization.h"
#include "math/bcl_math_mutate_result.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> MutateSheetOrder::s_Instance
    (
      GetObjectInstances().AddInstance( new MutateSheetOrder())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MutateSheetOrder::MutateSheetOrder() :
      m_ParallelProbability( 0.0)
    {
    }

    //! @brief constructor from a parallel conformation probability
    //! @param PARALLEL_PROBABILITY probability of ordering strands in a parallel conformation
    MutateSheetOrder::MutateSheetOrder( const double PARALLEL_PROBABILITY) :
      m_ParallelProbability( PARALLEL_PROBABILITY)
    {
      BCL_Assert
      (
        PARALLEL_PROBABILITY >= 0.0 && PARALLEL_PROBABILITY <= 1.0,
        "The given probability for ordering in parallel conformation should be between 0 and 1"
      );
    }

    //! @brief Clone function
    //! @return pointer to new MutateSheetOrder
    MutateSheetOrder *MutateSheetOrder::Clone() const
    {
      return new MutateSheetOrder( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &MutateSheetOrder::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &MutateSheetOrder::GetAlias() const
    {
      static const std::string s_name( "MutateSheetOrder");
      return s_name;
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer MutateSheetOrder::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "Orders strands in a sheet by sequence order.");
      serializer.AddInitializer
      (
        "parallel probability",
        "probability of ordering into parallel conformation",
        io::Serialization::GetAgent( &m_ParallelProbability)
      );

      return serializer;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief operator to sort the strands in the given sheet by sequence order
    //! @param SHEET Sheet to be mutated
    //! @return MutateResult that contains the mutated sheet
    math::MutateResult< assemble::Domain> MutateSheetOrder::operator()( const assemble::Domain &SHEET) const
    {
      // construct undefined sh ptr
      util::ShPtr< assemble::Domain> sp_undefined_domain;

      // make sure given domain is a sheet
      BCL_Assert
      (
        SHEET.GetTopology().IsDefined() &&
        (
          SHEET.GetTopology()->GetType() == assemble::Topology::e_Sheet ||
          SHEET.GetTopology()->GetType() == assemble::Topology::e_BetaBarrel
        ),
        "The given domain is not a sheet or a barrel"
      );

      // if there are less than two strands
      if( SHEET.GetNumberSSEs() < 2)
      {
        BCL_MessageVrb( "The given sheet has less than 2 strands");
        return math::MutateResult< assemble::Domain>( sp_undefined_domain, *this);
      }

      // determine to order parallel or not
      const bool order_parallel( random::GetGlobalRandom().Double() < m_ParallelProbability);

      // construct new sheet
      util::ShPtr< assemble::Domain> sp_sheet( new assemble::Domain());

      // get the subset of strands to be sorted and cast them to SSEs
      util::SiPtrVector< const assemble::SSE> strands( SHEET.GetTopology()->GetElements());

      // store the y flip
      // this is chosen instead of x flip to preserve the side chain orientations
      const coord::MoveRotateDefined y_flip( coord::MoveRotateDefined::GetFlipMove( coord::GetAxes().e_Y));

      // copy the first strands
      util::ShPtr< assemble::SSE> sp_first_strand( strands.FirstElement()->Clone());

      // decide to flip this or not randomly
      const bool flip_first( random::GetGlobalRandom().Boolean());

      // boolean whether previous strand was flipped or not
      bool is_previous_flipped( false);

      // iterate over strands
      for
      (
        util::SiPtrVector< const assemble::SSE>::const_iterator
          strand_itr( strands.Begin()), strand_itr_end( strands.End());
        strand_itr != strand_itr_end; ++strand_itr
      )
      {
        // boolean whether this strand should be flipped
        bool flip_this_strand( false);

        // if its the first SSE
        if( strand_itr == strands.Begin())
        {
          // then flip this only if flip_first was set to true
          flip_this_strand = flip_first;
        }
        else
        {
          // boolean whether this is parallel to the first one
          bool is_parallel
          (
            SHEET.GetTopology()->GetPackingForSSEGeometryPair
            (
              util::SiPtr< const assemble::SSEGeometryInterface>( *strand_itr),
              util::SiPtr< const assemble::SSEGeometryInterface>( *( strand_itr - 1))
            ).GetOrientation() ==
            assemble::SSEGeometryPacking::e_Parallel
          );

          // boolean whether this is parallel to previous one after changes to previous
          // it's still parallel only if it was parallel or the previous was flipped, but not both
          const bool is_parallel_after_changes( is_previous_flipped ^ is_parallel);

          // we want to flip if we're ordering anti-parallel and it's antiparallel to previous but not both
          flip_this_strand = order_parallel ^ is_parallel_after_changes;
        }

        // make a copy of the strand
        util::ShPtr< assemble::SSE> sp_new_strand( ( *strand_itr)->Clone());

        // if flip is requested
        if( flip_this_strand)
        {
          // apply flip
          y_flip.Move( *sp_new_strand);
        }

        // insert into sheet
        sp_sheet->Insert( sp_new_strand);

        // update boolean
        is_previous_flipped = flip_this_strand;
      }

      // sort and return the result
      return math::MutateResult< assemble::Domain>( sp_sheet, *this);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &MutateSheetOrder::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_ParallelProbability, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &MutateSheetOrder::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_ParallelProbability, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  } // namespace fold
} // namespace bcl

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
#include "fold/bcl_fold_mutate_sheet_sort.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_domain.h"
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
    const util::SiPtr< const util::ObjectInterface> MutateSheetSort::s_Instance
    (
      GetObjectInstances().AddInstance( new MutateSheetSort())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MutateSheetSort::MutateSheetSort() :
      m_ReverseSortProbability( 0.0)
    {
    }

    //! @brief constructor from a reverse sort probability
    //! @param REVERSE_SORT_PROBABILITY probability of sorting the strands in reverse order
    MutateSheetSort::MutateSheetSort
    (
      const double REVERSE_SORT_PROBABILITY
    ) :
      m_ReverseSortProbability( REVERSE_SORT_PROBABILITY)
    {
      BCL_Assert
      (
        REVERSE_SORT_PROBABILITY >= 0.0 && REVERSE_SORT_PROBABILITY <= 1.0,
        "The given probability for reverse sorting should be between 0 and 1"
      );
    }

    //! @brief Clone function
    //! @return pointer to new MutateSheetSort
    MutateSheetSort *MutateSheetSort::Clone() const
    {
      return new MutateSheetSort( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &MutateSheetSort::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &MutateSheetSort::GetAlias() const
    {
      static const std::string s_name( "MutateSheetSort");
      return s_name;
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer MutateSheetSort::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "Sorts strands in a sheet by sequence order.");
      serializer.AddInitializer
      (
        "reverse sort probability",
        "probability of sorting strands in reverse order",
        io::Serialization::GetAgent( &m_ReverseSortProbability)
      );

      return serializer;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief operator to sort the strands in the given sheet by sequence order
    //! @param SHEET Sheet to be mutated
    //! @return MutateResult that contains the mutated sheet
    math::MutateResult< assemble::Domain> MutateSheetSort::operator()( const assemble::Domain &SHEET) const
    {
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

      // determine to reverse order
      const bool reverse_order( random::GetGlobalRandom().Double() < m_ReverseSortProbability);

      // construct new sheet
      util::ShPtr< assemble::Domain> sp_sheet( new assemble::Domain());

      // get the subset of strands to be sorted and cast them to SSEs
      util::SiPtrVector< const assemble::SSE> strands( SHEET.GetTopology()->GetElements());
      // check that the cast worked
      BCL_Assert( strands.IsDefined(), "The dynamic cast failed from geometry to SSEs");

      // get the list of strands from domain which should be sorted by sequence order
      util::SiPtrVector< const assemble::SSE> strands_sorted( SHEET.GetSSEs());

      // if reverse order is requested
      if( reverse_order)
      {
        // reverse the sorted strands
        std::reverse( strands_sorted.Begin(), strands_sorted.End());
      }

      // construct a topology and update the elements for the new sheet
      util::ShPtr< assemble::Topology> sp_topology( new assemble::Topology());
      sp_topology->SetType( SHEET.GetTopology()->GetType());
      sp_topology->SetElements( strands_sorted);
      sp_topology->SetOrientationFromType();
      sp_sheet->SetTopology( sp_topology);

      // iterate strands and sorted strands at the same time
      for( size_t i( 0); i < strands.GetSize(); ++i)
      {
        // build transformation to move it to origin
        math::TransformationMatrix3D transformation( strands_sorted( i)->GetOrientation());
        transformation.Invert();

        // update the transformation to move it back to the sheet
        transformation( strands( i)->GetOrientation());

        // make a copy of the SSE
        util::ShPtr< assemble::SSE> sp_new_sse( strands_sorted( i)->Clone());

        // apply the transformation and add the transformed strand to sheet
        sp_new_sse->Transform( transformation);
        sp_sheet->Insert( sp_new_sse);
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
    std::istream &MutateSheetSort::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_ReverseSortProbability, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &MutateSheetSort::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_ReverseSortProbability, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  } // namespace fold
} // namespace bcl

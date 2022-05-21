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
#include "pdb/bcl_pdb_printer_score.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_quality_batch.h"
#include "pdb/bcl_pdb_line.h"
#include "score/bcl_score_protein_model_score_sum.h"
#include "util/bcl_util_sh_ptr_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace pdb
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> PrinterScore::s_Instance
    (
      GetObjectInstances().AddInstance( new PrinterScore())
    );

    //! @brief get format for table column format
    //! @return format for table column format
    const util::Format &PrinterScore::GetTableColumnFormat()
    {
      // initialize static const format
      static const util::Format s_format( util::Format().W( 13).FFP( 4).R());

      // end
      return s_format;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    PrinterScore::PrinterScore() :
      m_ScoringFunction(),
      m_Qualities()
    {
    }

    //! @brief construct with member variables
    //! @param SP_SCORE optional scoring function
    //! @param QUALITY_MEASURES set of quality measures to calculate
    PrinterScore::PrinterScore
    (
      const util::ShPtr< score::ProteinModelScoreSum> &SP_SCORE,
      const storage::Set< quality::Measure> &QUALITY_MEASURES
    ) :
      m_ScoringFunction( SP_SCORE),
      m_Qualities( QUALITY_MEASURES)
    {
    }

    //! @brief Clone function
    //! @return pointer to new PrinterScore
    PrinterScore *PrinterScore::Clone() const
    {
      return new PrinterScore( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &PrinterScore::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief writes a storage table to PDB lines
    //! @param TABLE table to be written
    //! @param WRITE_EMPTY_REMARK_LINE bool whether to precede table with empty remark line
    //! @return PDB lines
    util::ShPtrList< Line> PrinterScore::WriteTableToLines
    (
      const storage::Table< double> &TABLE,
      const bool WRITE_EMPTY_REMARK_LINE
    )
    {
      // initialize lines
      util::ShPtrList< Line> lines;

      // if the table is not empty
      if( !TABLE.IsEmpty())
      {
        // create an empty remark line
        util::ShPtr< Line> empty_remark_line( new Line( GetLineTypes().REMARK));
        empty_remark_line->Put( GetEntryTypes().REMARK_Number, util::Format()( size_t( s_RemarkNumber)));
        if( WRITE_EMPTY_REMARK_LINE)
        {
          lines.PushBack( empty_remark_line);
        }

        // write to stream
        std::stringstream ss;
        TABLE.WriteFormatted( ss, GetTableColumnFormat(), "Scores");
        storage::Vector< std::string> string_lines( util::StringLineListFromIStream( ss));

        // iterate through the lines
        for
        (
          storage::Vector< std::string>::const_iterator line_itr( string_lines.Begin()),
            line_itr_end( string_lines.End());
          line_itr != line_itr_end; ++line_itr
        )
        {
          // get a new line
          util::ShPtr< Line> new_line( empty_remark_line->Clone());

          // add the entry
          new_line->Put( GetEntryTypes().REMARK_String, *line_itr);
          lines.PushBack( new_line);
        }
      }

      // end
      return lines;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief takes a protein model and returns PDB lines
    //! @param PROTEIN_MODEL protein model to print
    //! @return PDB lines
    util::ShPtrList< Line> PrinterScore::operator()( const assemble::ProteinModel &PROTEIN_MODEL) const
    {
      // initialize table
      storage::Table< double> table
      (
        math::SumFunctionMixin< score::ProteinModel>::GetValueTableVerticalColumnNames()
      );

      // write score to file
      if( m_ScoringFunction.IsDefined())
      {
        table.Append( m_ScoringFunction->CreateValueTableVertical( PROTEIN_MODEL));
      }

      // add quality measures to table
      const assemble::QualityBatch quality_batch( m_Qualities, biol::GetAtomTypes().CA);
      if( !m_Qualities.IsEmpty()) // if qualities were given, add all the quality measures
      {
        table.Append( quality_batch.ConstructTable( PROTEIN_MODEL));
      }
      else // if no qualities were given, just add the completeness estimate
      {
        table.Append( quality_batch.ContructTableWithCompletenessOnly( PROTEIN_MODEL));
      }

      // return created lines
      return WriteTableToLines( table);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &PrinterScore::Read( std::istream &ISTREAM)
    {
      // read members

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &PrinterScore::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members

      // return the stream
      return OSTREAM;
    }

  } // namespace pdb
} // namespace bcl

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
#include "pdb/bcl_pdb_printer_loop_closure.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_locator_aa.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "fold/bcl_fold_locator_loop_domain.h"
#include "pdb/bcl_pdb_line.h"
#include "storage/bcl_storage_table.hpp"
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
    const util::SiPtr< const util::ObjectInterface> PrinterLoopClosure::s_Instance
    (
      GetObjectInstances().AddInstance( new PrinterLoopClosure())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    PrinterLoopClosure::PrinterLoopClosure() :
      m_ClosureThreshold( util::GetUndefined< double>())
    {
    }

    //! @brief constructor taking member variable parameters
    //! @param CLOSURE_THRESHOLD the threshold for considering the loop closed
    PrinterLoopClosure::PrinterLoopClosure( const double CLOSURE_THRESHOLD) :
      m_ClosureThreshold( CLOSURE_THRESHOLD)
    {
    }

    //! @brief Clone function
    //! @return pointer to new PrinterLoopClosure
    PrinterLoopClosure *PrinterLoopClosure::Clone() const
    {
      return new PrinterLoopClosure( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &PrinterLoopClosure::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief takes a protein model and returns PDB lines
    //! @param PROTEIN_MODEL protein model to print
    //! @return PDB lines
    util::ShPtrList< Line> PrinterLoopClosure::operator()( const assemble::ProteinModel &PROTEIN_MODEL) const
    {
      // initialize lines
      util::ShPtrList< Line> lines;

      // if no threshold was set
      if( !util::IsDefined( m_ClosureThreshold))
      {
        return lines;
      }

      // get the loop domain locator from model data
      util::ShPtr< util::ShPtrList< fold::LocatorLoopDomain> > sp_loop_domains
      (
        PROTEIN_MODEL.GetProteinModelData()->GetData( assemble::ProteinModelData::e_LoopDomainLocators)
      );

      // make sure it's defined
      if( sp_loop_domains.IsDefined())
      {
        util::Format format_num, format_short, format_med, format_long;
        format_num.W( 6).FFP( 3);
        format_short.W( 1);
        format_med.W( 4);
        format_long.W( 15);

        // create table which will hold all the information
        storage::Table< std::string> loop_info
        (
          storage::TableHeader
          (
            storage::Vector< std::string>::Create
            (
              format_short( "T"),
              format_long( "domain"),
              format_num( "length"),
              format_num( "omega"),
              format_num( "closed"),
              format_num( "rms"),
              format_med( "rms<" + util::Format()( m_ClosureThreshold))
            )
          )
        );

        // iterate through the loop domains to fill table with their loop closure status
        for
        (
          util::ShPtrList< fold::LocatorLoopDomain>::const_iterator
            loop_itr( sp_loop_domains->Begin()), loop_itr_end( sp_loop_domains->End());
          loop_itr != loop_itr_end;
          ++loop_itr
        )
        {
          const fold::LocatorLoopDomain &current_locator( **loop_itr);

          storage::Vector< std::string> row;
          assemble::LocatorAA left_aa_loc;
          assemble::LocatorAA right_aa_loc;
          if( current_locator.IsNToCSequenceDirection())
          {
            // amino acids in peptide bond
            left_aa_loc = current_locator.GetLoopSegments().LastElement().GetLocatorSSE().EndAALocator();
            right_aa_loc = assemble::LocatorAA( left_aa_loc.GetLocatorChain().GetChainID(), left_aa_loc.GetAAID() + 1, left_aa_loc.GetUsePDBID());

            row.PushBack( format_short( "C"));
          }
          else
          {
            // amino acids in peptide bond
            right_aa_loc = current_locator.GetLoopSegments().FirstElement().GetLocatorSSE().StartAALocator();
            left_aa_loc = assemble::LocatorAA( right_aa_loc.GetLocatorChain().GetChainID(), right_aa_loc.GetAAID() - 1, right_aa_loc.GetUsePDBID());

            row.PushBack( format_short( "N"));
          }

          const std::string flex_ident
          (
            util::Format()( current_locator.GetLoopSegments().FirstElement().GetLocatorSSE().GetSSEID().First()) +
            "-" +
            util::Format()( current_locator.GetLoopSegments().FirstElement().GetLocatorSSE().GetSSEID().Second()) +
            "," +
            util::Format()( current_locator.GetLoopSegments().LastElement().GetLocatorSSE().GetSSEID().First()) +
            "-" +
            util::Format()( current_locator.GetLoopSegments().LastElement().GetLocatorSSE().GetSSEID().Second())
          );
          row.PushBack( format_long( flex_ident));

          // pointer to left and right amino acid
          const util::SiPtr< const biol::AABase> aa_left( left_aa_loc.Locate( PROTEIN_MODEL));
          const util::SiPtr< const biol::AABase> aa_right( right_aa_loc.Locate( PROTEIN_MODEL));

          if( !aa_left.IsDefined() || !aa_right.IsDefined())
          {
            continue;
          }

          // row name using the amino acids
          const std::string row_name
          (
            util::Format()( aa_left->GetChainID()) + " " +
            util::Format()( aa_left->GetSeqID()) + " " + aa_left->GetType()->GetOneLetterCode() + " - " +
            util::Format()( aa_right->GetSeqID()) + " " + aa_right->GetType()->GetOneLetterCode()
          );

          // peptide bond parameters
          const storage::VectorND< 2, double> peptide_bond_length_angle( biol::AABase::PeptideBondLengthAndAngle( *aa_left, *aa_right));
          row.PushBack( format_num( peptide_bond_length_angle.First()));
          row.PushBack( format_num( peptide_bond_length_angle.Second()));
          row.PushBack( format_num( biol::AABase::AreAminoAcidsPeptideBonded( *aa_left, *aa_right, true)));

          // rms of pseudo residue
          const double rms( fold::LocatorLoopDomain::CalculateRMS( PROTEIN_MODEL, current_locator));
          row.PushBack( format_num( rms));
          row.PushBack( format_med( rms < m_ClosureThreshold));

          // add the row to the table of information
          loop_info.InsertRow( row_name, row, true);
        }

        // create an empty remark line
        util::ShPtr< Line> empty_remark_line( new Line( GetLineTypes().REMARK));
        empty_remark_line->Put( GetEntryTypes().REMARK_Number, util::Format()( size_t( s_RemarkNumber)));
        lines.PushBack( empty_remark_line);

        // write to stream
        std::stringstream ss;
        loop_info.WriteFormatted( ss, format_short, "Loops");
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
      else
      {
        BCL_MessageVrb( "Loop domain locators are not stored with the given model");
      }

      // convert table to lines
      return lines;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &PrinterLoopClosure::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_ClosureThreshold, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &PrinterLoopClosure::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_ClosureThreshold, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  } // namespace pdb
} // namespace bcl

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
#include "pdb/bcl_pdb_printer_membrane.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "biol/bcl_biol_membrane.h"
#include "pdb/bcl_pdb_line.h"
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
    const util::SiPtr< const util::ObjectInterface> PrinterMembrane::s_Instance
    (
      GetObjectInstances().AddInstance( new PrinterMembrane())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    PrinterMembrane::PrinterMembrane()
    {
    }

    //! @brief Clone function
    //! @return pointer to new PrinterMembrane
    PrinterMembrane *PrinterMembrane::Clone() const
    {
      return new PrinterMembrane( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &PrinterMembrane::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief takes a protein model and returns PDB lines
    //! @param PROTEIN_MODEL protein model to print
    //! @return PDB lines
    util::ShPtrList< Line> PrinterMembrane::operator()( const assemble::ProteinModel &PROTEIN_MODEL) const
    {
      // initialize lines
      util::ShPtrList< Line> lines;

      // get the membrane
      const util::ShPtr< biol::Membrane> sp_membrane
      (
        PROTEIN_MODEL.GetProteinModelData()->GetData( assemble::ProteinModelData::e_Membrane)
      );

      // if no membrane is found
      if( !sp_membrane.IsDefined())// && sp_membrane->IsDefined())
      {
        return lines;
      }

      // if membrane is undefined
      if( !sp_membrane->IsDefined())
      {
        return lines;
      }

      // create an empty remark line
      util::ShPtr< Line> empty_remark_line( new Line( GetLineTypes().REMARK));
      empty_remark_line->Put( GetEntryTypes().REMARK_Number, util::Format()( size_t( s_RemarkNumber)));
      lines.PushBack( empty_remark_line);

      // create format object
      util::Format number_format;
      number_format.FFP( 3);

      // write the normal
      util::ShPtr< Line> normal_line( empty_remark_line->Clone());
      normal_line->Put
      (
        GetEntryTypes().REMARK_String,
        "Membrane normal: " + number_format( sp_membrane->GetNormal().X()) + "\t" +
        number_format( sp_membrane->GetNormal().Y()) + "\t" + number_format( sp_membrane->GetNormal().Z())
      );
      lines.PushBack( normal_line);

      // write the center
      util::ShPtr< Line> center_line( empty_remark_line->Clone());
      center_line->Put
      (
        GetEntryTypes().REMARK_String,
        "Membrane center: " + number_format( sp_membrane->GetCenter().X()) + "\t" +
        number_format( sp_membrane->GetCenter().Y()) + "\t" + number_format( sp_membrane->GetCenter().Z())
      );
      lines.PushBack( center_line);

      // write the thickness
      util::ShPtr< Line> thickness_line( empty_remark_line->Clone());
      const linal::Vector< double> thicknesses( sp_membrane->GetThicknesses());
      std::string thickness_string( "Membrane thickness: ");
      for
      (
        linal::Vector< double>::const_iterator itr( thicknesses.Begin()), itr_end( thicknesses.End());
        itr != itr_end; ++itr
      )
      {
        thickness_string += number_format( *itr) + "\t";
      }
      thickness_line->Put( GetEntryTypes().REMARK_String, thickness_string);
      lines.PushBack( thickness_line);

      // end
      return lines;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &PrinterMembrane::Read( std::istream &ISTREAM)
    {
      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &PrinterMembrane::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // return the stream
      return OSTREAM;
    }

  } // namespace pdb
} // namespace bcl

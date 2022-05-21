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
#include "pdb/bcl_pdb_printer_biomatrix.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "assemble/bcl_assemble_protein_model_multiplier.h"
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
    const util::SiPtr< const util::ObjectInterface> PrinterBiomatrix::s_Instance
    (
      GetObjectInstances().AddInstance( new PrinterBiomatrix())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    PrinterBiomatrix::PrinterBiomatrix()
    {
    }

    //! @brief Clone function
    //! @return pointer to new PrinterBiomatrix
    PrinterBiomatrix *PrinterBiomatrix::Clone() const
    {
      return new PrinterBiomatrix( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &PrinterBiomatrix::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief takes a protein model and returns PDB lines
    //! @param PROTEIN_MODEL protein model to print
    //! @return PDB lines
    util::ShPtrList< Line> PrinterBiomatrix::operator()( const assemble::ProteinModel &PROTEIN_MODEL) const
    {
      // initialize lines to return
      util::ShPtrList< Line> bio_molecule_lines;

      // check that there is a multiplier
      const util::ShPtr< assemble::ProteinModelMultiplier> sp_multiplier
      (
        PROTEIN_MODEL.GetProteinModelData()->GetData( assemble::ProteinModelData::e_Multiplier)
      );
      if( !sp_multiplier.IsDefined())
      {
        // return empty list
        return bio_molecule_lines;
      }

      // create a map of ShPtr transformation matrices and chain ids
      storage::Map< char, storage::Vector< math::TransformationMatrix3D> > transformations;

      // iterate through the chain multipliers
      for
      (
        storage::Set< util::ShPtr< assemble::ChainMultiplier>, assemble::ChainMultiplierLessThan>::const_iterator
          multiplier_itr( sp_multiplier->GetChainMultipliers().Begin()),
          multiplier_itr_end( sp_multiplier->GetChainMultipliers().End());
        multiplier_itr != multiplier_itr_end; ++multiplier_itr
      )
      {
        // add the transformation to the map
        transformations[ ( *multiplier_itr)->GetInitialChainID()].PushBack
        (
          *( *multiplier_itr)->GetTransformationMatrix()
        );
      }

      // create an empty remark line
      util::ShPtr< Line> empty_remark_line( new Line( GetLineTypes().REMARK));
      empty_remark_line->Put( GetEntryTypes().REMARK_Number, util::Format()( size_t( s_RemarkNumber)));
      bio_molecule_lines.PushBack( empty_remark_line);

      // create biomolecule # line
      util::ShPtr< Line> biomolecule_nr_line( empty_remark_line->Clone());
      biomolecule_nr_line->Put( GetEntryTypes().REMARK_350_BiomoleculeIdentifier, "BIOMOLECULE: ");
      biomolecule_nr_line->Put( GetEntryTypes().REMARK_350_BiomoleculeNumber, "1");
      bio_molecule_lines.PushBack( biomolecule_nr_line);

      // iterate through the map of transformations
      for
      (
        storage::Map< char, storage::Vector< math::TransformationMatrix3D> >::const_iterator
          transformation_itr( transformations.Begin()), transformation_itr_end( transformations.End());
        transformation_itr != transformation_itr_end; ++transformation_itr
      )
      {
        // build the chain information line
        const std::string chain_text( "APPLY THE FOLLOWING TO CHAINS: ");
        util::ShPtr< Line> chain_line( empty_remark_line->Clone());
        chain_line->Put( GetEntryTypes().REMARK_String, chain_text + transformation_itr->first);
        bio_molecule_lines.PushBack( chain_line);

        // initialize transformation counter
        size_t transformation_ctr( 1);

        // write the transformation matrix to lines and append them
        for
        (
          storage::Vector< math::TransformationMatrix3D>::const_iterator
            matrix_itr( transformation_itr->second.Begin()), matrix_itr_end( transformation_itr->second.End());
          matrix_itr != matrix_itr_end; ++matrix_itr, ++transformation_ctr
        )
        {
          bio_molecule_lines.Append( WriteTransformationToLines( *matrix_itr, transformation_ctr));
        }
      }

      // end
      return bio_molecule_lines;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &PrinterBiomatrix::Read( std::istream &ISTREAM)
    {
      // read members

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &PrinterBiomatrix::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief write transformation matrix to Biomolecule lines
    //! @param TRANSFORMATION transformation matrix to be written
    //! @param COUNTER counter for current transformation
    //! @return formatted lines containing transformation information
     util::ShPtrList< Line> PrinterBiomatrix::WriteTransformationToLines
    (
      const math::TransformationMatrix3D &TRANSFORMATION,
      const size_t COUNTER
    )
    {
       // initialize lines to return
       util::ShPtrList< Line> bio_molecule_lines;

       // create an empty remark line
       util::ShPtr< Line> empty_remark_line( new Line( GetLineTypes().REMARK));
       empty_remark_line->Put( GetEntryTypes().REMARK_Number, "350");

       // get an iterator on the entry type
       EntryTypes::const_iterator type_itr( GetEntryTypes().REMARK_350_BIOMT1_Identifier.GetIterator());

       // iterate through the number of lines to add
       for( size_t i( 0); i != 3; ++i)
       {
         // create a new line with the correct identification
         util::ShPtr< Line> matrix_line_a( empty_remark_line->Clone());
         matrix_line_a->Put( *type_itr, "BIOMT" + util::Format()( i + 1));
         ++type_itr;
         matrix_line_a->Put( *type_itr, COUNTER);
         ++type_itr;

         // iterate through the rows of the matrix
         for( size_t j( 0); j != 4; ++j, ++type_itr)
         {
           // add the value to the line
           matrix_line_a->Put( *type_itr, TRANSFORMATION.GetMatrix()( j, i));
         }

         bio_molecule_lines.PushBack( matrix_line_a);
       }

       // end
       return bio_molecule_lines;
    }

  } // namespace pdb
} // namespace bcl

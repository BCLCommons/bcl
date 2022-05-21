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
#include "pdb/bcl_pdb_head.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_types.h"
#include "biol/bcl_biol_membrane.h"
#include "command/bcl_command_flag_interface.h"
#include "pdb/bcl_pdb_line_criterium.h"
#include "pdb/bcl_pdb_printer_membrane.h"
#include "pdb/bcl_pdb_residue_simple.h"
#include "storage/bcl_storage_triplet.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace pdb
  {

  ///////////
  // types //
  ///////////

    //! @brief operator for checking if one line is less than the other one
    //! @param LINE_LHS line on left-hand side
    //! @param LINE_RHS line on right-hand side
    //! @return if one line is less than the other one
    bool Head::RemarkLineLessThan::operator()
    (
      const util::ShPtr< Line> &LINE_LHS,
      const util::ShPtr< Line> &LINE_RHS
    ) const
    {
      // if both types are remark
      if
      (
        LINE_LHS.IsDefined() && LINE_RHS.IsDefined() &&
        LINE_LHS->GetType() == GetLineTypes().REMARK && LINE_RHS->GetType() == GetLineTypes().REMARK
      )
      {
        // compare remark #
        return LINE_LHS->GetString( GetEntryTypes().REMARK_Number) < LINE_RHS->GetString( GetEntryTypes().REMARK_Number);
      }

      // end
      return false;
    }

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> Head::s_Instance
    (
      GetObjectInstances().AddInstance( new Head())
    );

    //! BCL pdb identifier
    const std::string &Head::GetBCLPdbID()
    {
      static const std::string s_bcl_pdb_id( "BCL ");
      return s_bcl_pdb_id;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    Head::Head()
    {
    }

    //! @brief Clone function
    //! @return pointer to new Head
    Head *Head::Clone() const
    {
      return new Head( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Head::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief linetypes within group
    //! @return set of line types
    const storage::Set< LineType> &Head::GetTypesOfLines() const
    {
      static const storage::Set< LineType> s_line_types
      (
        GetLineTypes().HEADER.GetIterator(), GetLineTypes().MTRIX3.GetIterator() + 1
      );

      // end
      return s_line_types;
    }

    //! @brief access to lines of given type
    //! @param LINE_TYPE the desire line type
    //! @return lines of given type
    util::ShPtrList< Line> Head::GetLines( const LineType &LINE_TYPE) const
    {
      // search for that line type
      const storage::Map< LineType, util::ShPtrList< Line> >::const_iterator itr( m_Lines.Find( LINE_TYPE));

      if( itr != m_Lines.End())
      {
        return itr->second;
      }

      return util::ShPtrList< Line>();
    }

    //! @brief count the number of lines of given TYPE used for master record
    //! @param LINE_TYPE the line type
    //! @return the number of lines of that type found
    size_t Head::Count( const LineType &LINE_TYPE) const
    {
      // search for that line type
      const storage::Map< LineType, util::ShPtrList< Line> >::const_iterator itr( m_Lines.Find( LINE_TYPE));

      if( itr != m_Lines.End())
      {
        return itr->second.GetSize();
      }

      return 0;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief locate lines of given criterium
    //! @param CRITERIUM unary predicate with operator return true for line, that should be considered
    //! @return lines that are considered by criterium
    util::ShPtrList< Line> Head::CollectLines( const util::FunctionInterface< Line, bool> &CRITERIUM) const
    {
      util::ShPtrList< Line> lines;

      // iterate through map
      for
      (
        storage::Map< LineType, util::ShPtrList< Line> >::const_iterator itr( m_Lines.Begin()), itr_end( m_Lines.End());
        itr != itr_end;
        ++itr
      )
      {
        lines.Append( LineCriterium::Filter( itr->second, CRITERIUM));
      }

      // end
      return lines;
    }

    //! @brief locate lines of given criterium and line type
    //! @param CRITERIUM unary predicate with operator return true for line, that should be considered
    //! @param LINE_TYPE only consider this line type
    //! @return lines that are considered by criterium
    util::ShPtrList< Line>
    Head::CollectLines( const util::FunctionInterface< Line, bool> &CRITERIUM, const LineType &LINE_TYPE) const
    {
      // find lines of that type
      const storage::Map< LineType, util::ShPtrList< Line> >::const_iterator itr( m_Lines.Find( LINE_TYPE));

      // no lines of that type
      if( itr == m_Lines.End())
      {
        return util::ShPtrList< Line>();
      }

      // filter the lines according to the criterium
      return LineCriterium::Filter( itr->second, CRITERIUM);
    }

    //! @brief pushback a new line into that group
    //! @param ShPtr to the line
    //! @return true, if it fits into that group (line type is eligible)
    bool Head::PushBack( const util::ShPtr< Line> &LINE)
    {
      if( GetTypesOfLines().Contains( LINE->GetType()))
      {
        m_Lines[ LINE->GetType()].PushBack( LINE);
        return true;
      }

      return false;
    }

    //! @brief reset the line group
    void Head::Reset()
    {
      m_Lines.Reset();
    }

    //! extract all transformation matrices
    storage::Vector< math::TransformationMatrix3D> Head::GetTransformationMatrices() const
    {
      // get all the mtrix lines
      const util::ShPtrList< Line> &matrix_lines1( GetLines( GetLineTypes().MTRIX1));
      const util::ShPtrList< Line> &matrix_lines2( GetLines( GetLineTypes().MTRIX2));
      const util::ShPtrList< Line> &matrix_lines3( GetLines( GetLineTypes().MTRIX3));

      // holds the transformations
      storage::Vector< math::TransformationMatrix3D> transformationmatrices;

      // check that there is one rwo for matrix
      if( matrix_lines1.GetSize() != matrix_lines2.GetSize() || matrix_lines1.GetSize() != matrix_lines3.GetSize())
      {
        BCL_MessageCrt( "different number of MTRIXn lines!");
        return transformationmatrices;
      }

      // iterate over each row for all matrices
      for
      (
        util::ShPtrList< Line>::const_iterator
          matrix1_line_itr( matrix_lines1.Begin()), matrix1_line_itr_end( matrix_lines1.End()),
          matrix2_line_itr( matrix_lines2.Begin()), matrix2_line_itr_end( matrix_lines2.End()),
          matrix3_line_itr( matrix_lines3.Begin()), matrix3_line_itr_end( matrix_lines3.End());
        matrix1_line_itr != matrix1_line_itr_end &&
        matrix2_line_itr != matrix2_line_itr_end &&
        matrix3_line_itr != matrix3_line_itr_end;
        ++matrix1_line_itr, ++matrix2_line_itr, ++matrix3_line_itr
      )
      {
        // construct matrix
        linal::Matrix< double> new_matrix( 4, 4, double( 0));

        // first col
        new_matrix( 0, 0) = ( *matrix1_line_itr)->GetNumericalValue< double>( GetEntryTypes().MTRIX1_M11);
        new_matrix( 1, 0) = ( *matrix1_line_itr)->GetNumericalValue< double>( GetEntryTypes().MTRIX1_M12);
        new_matrix( 2, 0) = ( *matrix1_line_itr)->GetNumericalValue< double>( GetEntryTypes().MTRIX1_M13);
        new_matrix( 3, 0) = ( *matrix1_line_itr)->GetNumericalValue< double>( GetEntryTypes().MTRIX1_V1);

        // second col
        new_matrix( 0, 1) = ( *matrix2_line_itr)->GetNumericalValue< double>( GetEntryTypes().MTRIX2_M21);
        new_matrix( 1, 1) = ( *matrix2_line_itr)->GetNumericalValue< double>( GetEntryTypes().MTRIX2_M22);
        new_matrix( 2, 1) = ( *matrix2_line_itr)->GetNumericalValue< double>( GetEntryTypes().MTRIX2_M23);
        new_matrix( 3, 1) = ( *matrix2_line_itr)->GetNumericalValue< double>( GetEntryTypes().MTRIX2_V2);

        // last col
        new_matrix( 0, 2) = ( *matrix3_line_itr)->GetNumericalValue< double>( GetEntryTypes().MTRIX3_M31);
        new_matrix( 1, 2) = ( *matrix3_line_itr)->GetNumericalValue< double>( GetEntryTypes().MTRIX3_M32);
        new_matrix( 2, 2) = ( *matrix3_line_itr)->GetNumericalValue< double>( GetEntryTypes().MTRIX3_M33);
        new_matrix( 3, 2) = ( *matrix3_line_itr)->GetNumericalValue< double>( GetEntryTypes().MTRIX3_V3);

        // insert
        transformationmatrices.PushBack( math::TransformationMatrix3D( new_matrix));
      }

      // end
      return transformationmatrices;
    }

    //! read one column of a BIOMT matrix remark 350 line into a matrix
    //! @param LINE the pdb line to read from
    //! @param BIOMT_FIRST_COL_ENTRY_TYPE any of REMARK_350_BIOMT1_M11, REMARK_350_BIOMT1_M21, REMARK_350_BIOMT1_M31
    //! @param COL_NR which col is read
    //! @param MATRIX the matrix values are added to
    void ReadREMARK_350_BIOMTMatrixCol
    (
      const Line &LINE,
      const EntryType &BIOMT_FIRST_COL_ENTRY_TYPE,
      const size_t COL_NR,
      linal::Matrix< double> &MATRIX
    )
    {
      EntryTypes::const_iterator itr( BIOMT_FIRST_COL_ENTRY_TYPE.GetIterator());
      MATRIX( 0, COL_NR) = LINE.GetNumericalValue< double>( *itr);
      MATRIX( 1, COL_NR) = LINE.GetNumericalValue< double>( *++itr);
      MATRIX( 2, COL_NR) = LINE.GetNumericalValue< double>( *++itr);
      MATRIX( 3, COL_NR) = LINE.GetNumericalValue< double>( *++itr);
    }

    //! @brief the transformation matrices to generate different bio molecules by applying transformations to different chains
    //! @param CHAIN_IDS all chain ids that should be considered
    //! @return Map of biomolecule number to a vector of transformations of chainid it should be applied to, the new chain id and the transformation
    storage::Map< size_t, storage::Vector< storage::Triplet< char, char, math::TransformationMatrix3D> > >
    Head::GetBioTransformationMatrices( const std::string &CHAIN_IDS) const
    {
      storage::Map< size_t, storage::Vector< storage::Triplet< char, char, math::TransformationMatrix3D> > >
        transformation_matrices;

      // search for all REMARK  350 lines
      LineCriterium criterium;
      criterium.AddCriterium( GetEntryTypes().REMARK_Number, 350);

      // acquire matrix lines
      const util::ShPtrList< Line> matrix_lines( CollectLines( criterium, GetLineTypes().REMARK));

      const std::string biomolecule_identifier( "BIOMOLECULE:");
      const std::string chains_identifier( "CHAINS:");

      // if there are no chain ids
      if( CHAIN_IDS.empty())
      {
        return transformation_matrices;
      }
      char next_chain_id( *( --CHAIN_IDS.end()));

      // create set of transformed chains
      storage::Set< char> transformed_chains;

      util::ShPtrList< Line>::const_iterator line_itr( matrix_lines.Begin()), line_itr_end( matrix_lines.End());

      // iterate over all lines
      while( line_itr != line_itr_end)
      {
        // check for BIOMOLECULE string
        while
        (
             line_itr != line_itr_end
          && ( *line_itr)->GetString( GetEntryTypes().REMARK_350_BiomoleculeIdentifier) != biomolecule_identifier
        )
        {
          ++line_itr;
        }
        if( line_itr == line_itr_end)
        {
          BCL_MessageVrb( "could not find any more " + biomolecule_identifier);
          break;
        }

        // report that BIOMOLECULE was found
        const size_t current_biomolecule
        (
          ( *line_itr)->GetNumericalValue< size_t>( GetEntryTypes().REMARK_350_BiomoleculeNumber)
        );
        BCL_MessageStd( "found " + biomolecule_identifier + ' ' + util::Format()( current_biomolecule));

        // goto next line
        ++line_itr;

        // check for apply identifier
        while( line_itr != line_itr_end)
        {
          // break if new biomolecule appears
          if( ( *line_itr)->GetString( GetEntryTypes().REMARK_350_BiomoleculeIdentifier) == biomolecule_identifier)
          {
            break;
          }

          // check that the line contains a chains_identifier
          if( ( *line_itr)->GetString( GetEntryTypes().REMARK_350_ChainsIdentifier) != chains_identifier)
          {
            ++line_itr;
            continue;
          }

          // extract the chains the transformation has to be applied to
          std::string chains_string( ( *line_itr)->GetString( GetEntryTypes().REMARK_350_ChainsList));

          while( ++line_itr != line_itr_end && ( *line_itr)->GetString( GetEntryTypes().REMARK_350_ChainsIdentifier) == chains_identifier)
          {
            chains_string += ( *line_itr)->GetString( GetEntryTypes().REMARK_350_ChainsList);
          }

          // sort the chain ids and erase every non alpha character
          std::sort( chains_string.begin(), chains_string.end());
          chains_string.erase
          (
            std::remove_if( chains_string.begin(), chains_string.end(), std::not1( std::ptr_fun( &isalpha))),
            chains_string.end()
          );

          // remove all chain id, that are not in the argument
          {
            std::string arg_chainids( CHAIN_IDS);
            std::sort( arg_chainids.begin(), arg_chainids.end());
            std::string filtered_chainids;
            std::set_intersection
            (
              chains_string.begin(), chains_string.end(),
              arg_chainids.begin(), arg_chainids.end(),
              std::inserter( filtered_chainids, filtered_chainids.begin())
            );

            chains_string = filtered_chainids;
          }

          BCL_MessageCrt
          (
            "Chains to apply transformations to:\n" + util::Format()( chains_string)
          );

          // should be at BIOMT1
          // read all matrices (one matrix consists of BIMT1, 2 and 3)
          while
          (
               line_itr != line_itr_end
            && ( *line_itr)->GetString( GetEntryTypes().REMARK_350_BIOMT1_Identifier) == "BIOMT1"
          )
          {
            linal::Matrix< double> new_matrix( 4, 4, double( 0));
            new_matrix( 3, 3) = 1;

            BCL_Assert
            (
              ( *line_itr)->GetString( GetEntryTypes().REMARK_350_BIOMT1_Identifier) == "BIOMT1",
              "wrong identifier; should be \"BIOMT1\", but line is: " + ( *line_itr)->GetString()
            );
            ReadREMARK_350_BIOMTMatrixCol( **line_itr, GetEntryTypes().REMARK_350_BIOMT1_M11, 0, new_matrix);

            //go to next line
            ++line_itr;
            BCL_Assert( line_itr != line_itr_end, "ran out of lines");
            BCL_Assert
            (
              ( *line_itr)->GetString( GetEntryTypes().REMARK_350_BIOMT2_Identifier) == "BIOMT2",
              "wrong identifier; should be \"BIOMT2\", but line is: " + ( *line_itr)->GetString()
            );
            ReadREMARK_350_BIOMTMatrixCol( **line_itr, GetEntryTypes().REMARK_350_BIOMT2_M21, 1, new_matrix);

            //go to next line
            ++line_itr;
            BCL_Assert( line_itr != line_itr_end, "ran out of lines");
            BCL_Assert
            (
              ( *line_itr)->GetString( GetEntryTypes().REMARK_350_BIOMT3_Identifier) == "BIOMT3",
              "wrong identifier; should be \"BIOMT3\", but line is: " + ( *line_itr)->GetString()
            );
            ReadREMARK_350_BIOMTMatrixCol( **line_itr, GetEntryTypes().REMARK_350_BIOMT3_M31, 2, new_matrix);

            // iterate over all chains this matrix will be applied to
            for
            (
              std::string::const_iterator itr( chains_string.begin()), itr_end( chains_string.end());
              itr != itr_end; ++itr
            )
            {
              // initialize next chain id
              char new_chain_id;

              // if this chain has already been transformed
              if( transformed_chains.Contains( *itr))
              {
                // assign to next chain
                do
                {
                  ++next_chain_id;
                }
                while( !isalpha( next_chain_id));

                new_chain_id = next_chain_id;
              }
              // this is the first transformation for this chain
              else
              {
                // the new chain will have the same chain id
                new_chain_id = *itr;
                transformed_chains.Insert( *itr);
              }

              // insert the transformation matrix and the corresponding chains it should be applied
              transformation_matrices[ current_biomolecule].PushBack
              (
                storage::Triplet< char, char, math::TransformationMatrix3D>
                (
                  *itr,
                  new_chain_id,
                  math::TransformationMatrix3D( new_matrix)
                )
              );
            }
            //go to next line
            ++line_itr;
          } // loop over matrices
        } // loop over "apply to chain"
      } // loop over biomolecules

      return transformation_matrices;
    }

    //! @brief get the membrane from the remark lines
    //! @return the membrane from the remark lines
    util::ShPtr< biol::Membrane> Head::GetMembrane() const
    {
      // search for all membrane lines
      LineCriterium criterium;
      criterium.AddCriterium( GetEntryTypes().REMARK_Number, size_t( PrinterMembrane::s_RemarkNumber));

      // acquire matrix lines
      const util::ShPtrList< Line> membrane_lines( CollectLines( criterium, GetLineTypes().REMARK));

      // if no membrane lines or wrong size
      if( membrane_lines.GetSize() != 4)
      {
        return util::ShPtr< biol::Membrane>();
      }

      // skip the first line
      util::ShPtrList< Line>::const_iterator line_itr( membrane_lines.Begin());
      ++line_itr;

      // second line is normal
      const storage::Vector< std::string> normal_entries( util::SplitString( ( *line_itr)->GetString()));
      BCL_Assert( normal_entries.GetSize() == 7, "There should be 7 membrane normal lines in the PDB");
      const linal::Vector3D normal
      (
        util::ConvertStringToNumericalValue< double>( normal_entries( 4)),
        util::ConvertStringToNumericalValue< double>( normal_entries( 5)),
        util::ConvertStringToNumericalValue< double>( normal_entries( 6))
      );
      ++line_itr;

      // third line is center
      const storage::Vector< std::string> center_entries( util::SplitString( ( *line_itr)->GetString()));
      BCL_Assert( center_entries.GetSize() == 7, "There should be 7 membrane center lines in the PDB");
      const linal::Vector3D center
      (
        util::ConvertStringToNumericalValue< double>( center_entries( 4)),
        util::ConvertStringToNumericalValue< double>( center_entries( 5)),
        util::ConvertStringToNumericalValue< double>( center_entries( 6))
      );
      ++line_itr;

      // last line is thickness, use command line thickness if set
      storage::Vector< double> thickness;

      if( biol::Membrane::GetFlagMembrane()->GetFirstParameter()->GetWasSetInCommandLine())
      {
        thickness = biol::Membrane::GetCommandLineMembrane().GetThicknesses();
      }
      else
      {
        const storage::Vector< std::string> thickness_entries( util::SplitString( ( *line_itr)->GetString()));
        BCL_Assert( thickness_entries.GetSize() >= 9, "There should be 9 membrane thickness lines in the PDB");
        thickness = storage::Vector< double>::Create
        (
          util::ConvertStringToNumericalValue< double>( thickness_entries( 4)),
          util::ConvertStringToNumericalValue< double>( thickness_entries( 5)),
          util::ConvertStringToNumericalValue< double>( thickness_entries( 6)),
          util::ConvertStringToNumericalValue< double>( thickness_entries( 7)),
          util::ConvertStringToNumericalValue< double>( thickness_entries( 8))
        );
      }

      // end
      return util::ShPtr< biol::Membrane>( new biol::Membrane( thickness, normal, center));
    }

    //! @brief gets the pdb ID read in from the pdb file
    //! @return the pdb ID read in from the pdb file
    std::string Head::GetPDBID() const
    {
      // find header lines
      const util::ShPtrList< Line> header_lines( GetLines( GetLineTypes().HEADER));

      //
      if( header_lines.IsEmpty())
      {
        return GetBCLPdbID();
      }

      // get pdb id from header line
      const std::string pdb_id( header_lines.FirstElement()->GetString( GetEntryTypes().HEADERIDCode));

      // check that it is not empty
      if( util::TrimString( pdb_id).empty())
      {
        return GetBCLPdbID();
      }

      // return the pdb id
      return pdb_id;
    }

    //! @brief full name for all het residues
    //! @return map with key res name and data fullname
    storage::Map< std::string, std::string> Head::GetHetFullname() const
    {
      storage::Map< std::string, std::string> het_fullname;

      // locate all lines HETNAM
      const storage::Map< LineType, util::ShPtrList< Line> >::const_iterator find_itr( m_Lines.Find( GetLineTypes().HETNAM));
      if( find_itr == m_Lines.End())
      {
        return het_fullname;
      }

      // store the current line number for each HETNAM short name
      storage::Map< std::string, size_t> nam_lines;

      // iterate over lines
      for
      (
        util::ShPtrList< Line>::const_iterator itr( find_itr->second.Begin()), itr_end( find_itr->second.End());
        itr != itr_end;
        ++itr
      )
      {
        const std::string current( ( *itr)->GetString( GetEntryTypes().HETNAMIdentifier));
        const size_t line_count( ++nam_lines[ current]);
        if( line_count > 1)
        {
          const size_t current_linenumber( ( *itr)->GetNumericalValue< size_t>( GetEntryTypes().HETNAMContinuation));
          // if the current line number is not defined, it would be the second one without a line number, which is implicitly '1'
          // or if both line numbers are the same
          if( !util::IsDefined( current_linenumber) || current_linenumber != line_count)
          {
            BCL_MessageCrt
            (
              "found HETNAM line that does not fit in continuation count: " + ( *itr)->GetString()
            )
            continue;
          }
        }
        het_fullname[ current].append( util::TrimString( ( *itr)->GetString( GetEntryTypes().HETNAMText)));
      }

      // end
      return het_fullname;
    }

    //! @brief formula for all het residues
    //! @return map with key res name and data formula
    storage::Map< std::string, std::string> Head::GetHetFormula() const
    {
      storage::Map< std::string, std::string> het_formula;

      // locate all lines FORMUL
      const storage::Map< LineType, util::ShPtrList< Line> >::const_iterator find_itr( m_Lines.Find( GetLineTypes().FORMUL));
      if( find_itr == m_Lines.End())
      {
        return het_formula;
      }

      // iterate over lines
      for
      (
        util::ShPtrList< Line>::const_iterator itr( find_itr->second.Begin()), itr_end( find_itr->second.End());
        itr != itr_end;
        ++itr
      )
      {
        het_formula[ ( *itr)->GetString( GetEntryTypes().FORMULIdentifier)].append( ( *itr)->GetString( GetEntryTypes().FORMULChemicalFormula));
      }

      // end
      return het_formula;
    }

    //! @brief create the chains as they were given in the seqres
    //! @return map with chainid as key and list of residues as data
    storage::Map< char, storage::List< ResidueSimple> > Head::GetSEQRESProteinChains() const
    {
      storage::Map< char, storage::List< ResidueSimple> > chains;

      // retrieve seqres lines for all chains
      const util::ShPtrList< Line> seqreslines( GetLines( GetLineTypes().SEQRES));

      // store number of residues for that chain
      storage::Map< char, size_t> number_res_in_chain;

      // iterate over pdb seqres lines
      for
      (
        util::ShPtrList< Line>::const_iterator seq_line_itr( seqreslines.Begin()), seq_line_itr_end( seqreslines.End());
        seq_line_itr != seq_line_itr_end;
        ++seq_line_itr
      )
      {
        const char current_chain_id( ( *seq_line_itr)->GetChar( GetEntryTypes().SEQRESChainID));
        storage::Map< char, size_t>::const_iterator num_res_itr( number_res_in_chain.Find( current_chain_id));
        // check if there is already a number of residues stored for that chain
        if( num_res_itr == number_res_in_chain.End())
        {
          // store number of residues for that chain
          number_res_in_chain[ current_chain_id] =
            ( *seq_line_itr)->GetNumericalValue< size_t>( GetEntryTypes().SEQRESNrOfResiduesInChain);
        }
        else
        {
          // assert that all seqres lines for that chain have the same number of residues
          BCL_Assert
          (
            num_res_itr->second ==
              ( *seq_line_itr)->GetNumericalValue< size_t>( GetEntryTypes().SEQRESNrOfResiduesInChain),
            "number of residues in seqreslines for that chain are not all equal: \'" + util::Format()( current_chain_id)
            + "\'"
          );
        }

        // iterate over all 13 residues in a seqres line
        for
        (
          EntryTypes::const_iterator
            res_itr( GetEntryTypes().SEQRESName_1.GetIterator()),
            res_itr_end( ++GetEntryTypes().SEQRESName_13.GetIterator());
          res_itr != res_itr_end;
          ++res_itr
        )
        {
          const std::string res_name( util::TrimString( ( *seq_line_itr)->GetString( *res_itr)));

          // is this residue entry given in this SEQRES line
          if( res_name.size() < 3)
          {
            break;
          }

          chains[ current_chain_id].PushBack( ResidueSimple( res_name, current_chain_id));
        }
      }

      // check that number of res in seqresline SequenceNrOfResiduesInChain does match the actual number of residues
      for
      (
        storage::Map< char, size_t>::const_iterator
          itr( number_res_in_chain.Begin()), itr_end( number_res_in_chain.End());
        itr != itr_end;
        ++itr
      )
      {
        // find the chain with residues for that chain id
        const storage::Map< char, storage::List< ResidueSimple> >::const_iterator chain_itr( chains.Find( itr->first));

        // check if there have been residues for that chain at all
        if( chain_itr == chains.End())
        {
          BCL_MessageCrt
          (
            "Although SEQRES lines for this chain were found: \'" + util::Format()( itr->first) + "\' no residues " +
            "have been recognized. Probably it is a DNA or RNA chain which is not handled yet"
          );
        }

        // check that number of residues for that chain match
        else if( chain_itr->second.GetSize() != itr->second)
        {
          BCL_MessageCrt
          (
            "number of pretended residues in SEQRES lines for that chain \'" + util::Format()( itr->first) + "\' " +
            "does not match with actual number of residues: " + util::Format()( itr->second) + " != " +
            util::Format()( chains[ itr->first].GetSize())
          );
        }
      }

      // end
      return chains;
    }

    //! @brief map of residues that could not be located in the experiments using REMARK 465
    //! @param MODEL_NUMBER missing residues for the given model number
    //! @return map of chains and ShPtrList of Residues that are not in the ATOM lines
    storage::Map< char, storage::List< ResidueSimple> > Head::GetMissingResidues( const size_t MODEL_NUMBER) const
    {
      // residues
      storage::Map< char, storage::List< ResidueSimple> > residues;

      // search for all REMARK  465 lines
      LineCriterium criteria;
      criteria.SetMeetAllCriteria( true);
      criteria.AddCriterium( GetEntryTypes().REMARK_Number, 465);

      // acquire missing residue lines
      const util::ShPtrList< Line> remark_465_lines( CollectLines( criteria, GetLineTypes().REMARK));
      util::ShPtrList< Line> remark_465_lines_no_header;
      {
        // header identifier
        const std::string header_identifier( "M RES C SSSEQI");
        util::ShPtrList< Line>::const_iterator line_itr( remark_465_lines.Begin()), line_itr_end( remark_465_lines.End());
        // check for header string
        while
        (
             line_itr != line_itr_end
          && util::TrimString( ( *line_itr)->GetString( GetEntryTypes().REMARK_String)) != header_identifier
        )
        {
          ++line_itr;
        }
        if( line_itr != line_itr_end)
        {
          ++line_itr;
          remark_465_lines_no_header = util::ShPtrList< Line>( line_itr, line_itr_end);
        }
      }

      criteria.AddCriterium( GetEntryTypes().REMARK_465_ModelNumber, MODEL_NUMBER);
      util::ShPtrList< Line> missing_res_lines( LineCriterium::Filter( remark_465_lines_no_header, criteria));

      if( missing_res_lines.IsEmpty() && MODEL_NUMBER == 1)
      {
        criteria.Reset();
        criteria.AddCriterium( GetEntryTypes().REMARK_465_ModelNumber, "");
        missing_res_lines = LineCriterium::Filter( remark_465_lines_no_header, criteria);
      }

      if( missing_res_lines.IsEmpty())
      {
        return residues;
      }

      // as long as their are remaining lines
      for
      (
        util::ShPtrList< Line>::const_iterator
          line_itr( missing_res_lines.Begin()), line_itr_end( missing_res_lines.End());
        line_itr != line_itr_end;
        ++line_itr
      )
      {
        // get residue name
        const std::string residue_name( ( *line_itr)->GetString( GetEntryTypes().REMARK_465_ResidueName));

        // check that it is valid residue name
        if( !biol::GetAATypes().AATypeFromThreeLetterCode( residue_name).IsDefined())
        {
          BCL_MessageCrt( "found undefined residues in missing residue line: " + ( *line_itr)->GetString());
          // end
          return residues;
        }

        const char chain_id( ( *line_itr)->GetChar( GetEntryTypes().REMARK_465_ChainID));
        // current residue
        const ResidueSimple current_residue
        (
          residue_name,
          chain_id,
          ( *line_itr)->GetNumericalValue< int>( GetEntryTypes().REMARK_465_ResidueSequenceID),
          ( *line_itr)->GetChar( GetEntryTypes().REMARK_465_InsertionCode)
        );

        // insert
        residues[ chain_id].PushBack( current_residue);
      }

      // end
      return residues;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @return outputstream which was written to
    std::ostream &Head::WriteLines( std::ostream &OSTREAM) const
    {
      // iterate through map
      for
      (
        storage::Map< LineType, util::ShPtrList< Line> >::const_iterator itr( m_Lines.Begin()), itr_end( m_Lines.End());
        itr != itr_end;
        ++itr
      )
      {
        // if this is a remark section
        if( itr->first == GetLineTypes().REMARK)
        {
          // copy the lines
          util::ShPtrList< Line> remark_lines( itr->second);

          // sort by remark #
          remark_lines.Sort( RemarkLineLessThan());

          // iterate through all lines
          for
          (
            util::ShPtrList< Line>::const_iterator line_itr( remark_lines.Begin()), line_itr_end( remark_lines.End());
            line_itr != line_itr_end;
            ++line_itr
          )
          {
            OSTREAM << ( *line_itr)->GetString() << '\n';
          }
        }
        else
        {
          // iterate through all lines
          for
          (
            util::ShPtrList< Line>::const_iterator line_itr( itr->second.Begin()), line_itr_end( itr->second.End());
            line_itr != line_itr_end;
            ++line_itr
          )
          {
            OSTREAM << ( *line_itr)->GetString() << '\n';
          }
        }
      }

      // end
      return OSTREAM;
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &Head::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Lines, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &Head::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Lines, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief pdb line to sse definition form a HELIX or SHEET line, with sstype, start and end residue
    //! @return Triplet of SSType starting and end residue
    storage::Triplet< biol::SSType, ResidueSimple, ResidueSimple>
    Head::SSEDefinitionFromPDBLine( const Line &PDB_LINE)
    {
      // helix definition
      if( PDB_LINE.GetType() == GetLineTypes().HELIX)
      {
        return
          storage::Triplet< biol::SSType, ResidueSimple, ResidueSimple>
          (
            biol::GetSSTypes().SSTypeFromPDBHelixClass( PDB_LINE.GetNumericalValue< size_t>( GetEntryTypes().HELIXClass)),
            ResidueSimple
            (
              PDB_LINE.GetString( GetEntryTypes().HELIXResidueName_Initial),
              PDB_LINE.GetChar( GetEntryTypes().HELIXChainID_Initial),
              PDB_LINE.GetNumericalValue< int>( GetEntryTypes().HELIXSequenceID_Initial),
              PDB_LINE.GetChar( GetEntryTypes().HELIXInsertionCode_Initial)
            ),
            ResidueSimple
            (
              PDB_LINE.GetString( GetEntryTypes().HELIXResidueName_Terminal),
              PDB_LINE.GetChar( GetEntryTypes().HELIXChainID_Terminal),
              PDB_LINE.GetNumericalValue< int>( GetEntryTypes().HELIXSequenceID_Terminal),
              PDB_LINE.GetChar( GetEntryTypes().HELIXInsertionCode_Terminal)
            )
          );
      }
      // sheet definition
      else if( PDB_LINE.GetType() == GetLineTypes().SHEET)
      {
        return
          storage::Triplet< biol::SSType, ResidueSimple, ResidueSimple>
          (
            biol::GetSSTypes().STRAND,
            ResidueSimple
            (
              PDB_LINE.GetString( GetEntryTypes().SHEETResidueName_Initial),
              PDB_LINE.GetChar( GetEntryTypes().SHEETChainID_Initial),
              PDB_LINE.GetNumericalValue< int>( GetEntryTypes().SHEETSequenceID_Initial),
              PDB_LINE.GetChar( GetEntryTypes().SHEETInsertionCode_Initial)
            ),
            ResidueSimple
            (
              PDB_LINE.GetString( GetEntryTypes().SHEETResidueName_Terminal),
              PDB_LINE.GetChar( GetEntryTypes().SHEETChainID_Terminal),
              PDB_LINE.GetNumericalValue< int>( GetEntryTypes().SHEETSequenceID_Terminal),
              PDB_LINE.GetChar( GetEntryTypes().SHEETInsertionCode_Terminal)
            )
          );
      }

      // cannot handle other sse definitions yet
      BCL_MessageCrt( "cannot generate sse definition for that pdbline\n" + PDB_LINE.GetString());

      // end
      return storage::Triplet< biol::SSType, ResidueSimple, ResidueSimple>();
    }

    //! @brief get all sse definitions as they are given in the HELIX/SHEET section
    //! @param HELIX_CLASSES set of helix types to consider
    //! @param MERGE_OVERLAPPING merge overlapping sses of the same type
    //! @return Map of chain ids and
    storage::Map< char, storage::List< storage::Triplet< biol::SSType, ResidueSimple, ResidueSimple> > >
    Head::SSEDefinitions( const storage::Set< biol::SSType> &HELIX_CLASSES, const bool MERGE_OVERLAPPING) const
    {
      storage::Map< char, storage::List< storage::Triplet< biol::SSType, ResidueSimple, ResidueSimple> > >
        set_of_sse_definitions;

      // extract HELIX and SHEET lines form pdb
      util::ShPtrList< Line> helix_strand_lines;
      helix_strand_lines.Append( GetLines( GetLineTypes().HELIX));
      helix_strand_lines.Append( GetLines( GetLineTypes().SHEET));

      // insert helix and sheet information from pdb into ordered set
      for
      (
        util::ShPtrList< Line>::const_iterator itr( helix_strand_lines.Begin()), itr_end( helix_strand_lines.End());
        itr != itr_end; ++itr
      )
      {
        // create the sse definition
        const storage::Triplet< biol::SSType, ResidueSimple, ResidueSimple>
          sse_definition( SSEDefinitionFromPDBLine( **itr));

        const char chain_id( sse_definition.Second().GetChainID());
        // check if starting and terminating residue are from the same chain
        if( chain_id != sse_definition.Third().GetChainID())
        {
          BCL_MessageCrt( "start and end residue for that line are from different chains\n" + ( *itr)->GetString());
          continue; // next sse definition
        };

        // check that residues are given in correct order
        if( sse_definition.Third() < sse_definition.Second())
        {
          BCL_MessageCrt
          (
            "SSE definition in this line might be given in wrong order, so skipping!\n" + ( *itr)->GetString()
          );
          continue;
        }

        // check that helix types are considered
        if( ( *itr)->GetType() == GetLineTypes().HELIX && HELIX_CLASSES.Find( sse_definition.First()) == HELIX_CLASSES.End())
        {
          BCL_MessageStd( "HELIX line is ignored, since helix class is not considered:\n" + ( *itr)->GetString());
          continue;
        }

        set_of_sse_definitions[ chain_id].PushBack( sse_definition);
      }

      if( MERGE_OVERLAPPING)
      {
        return SSEDefinitionsMergeOverlap( set_of_sse_definitions);
      }

      return set_of_sse_definitions;
    }

    //! @brief check and merge sse definitions for overlap
    //! @param SSE_DEFINITIONS
    //! @return Map of chain ids and sse definitions that were merged when possible
    storage::Map< char, storage::List< storage::Triplet< biol::SSType, ResidueSimple, ResidueSimple> > >
    Head::SSEDefinitionsMergeOverlap( const storage::Map< char, storage::List< storage::Triplet< biol::SSType, ResidueSimple, ResidueSimple> > > &SSE_DEFINITIONS)
    {
      typedef std::pointer_to_binary_function< const storage::Triplet< biol::SSType, ResidueSimple, ResidueSimple> &, const storage::Triplet< biol::SSType, ResidueSimple, ResidueSimple> &, bool> SSEInfoCompareType;

      // initiate set of sse information
      std::map< char, std::set< storage::Triplet< biol::SSType, ResidueSimple, ResidueSimple>, SSEInfoCompareType> >
        set_of_sse_definitions;

      // insert helix and sheet information from pdb into ordered set
      for
      (
        storage::Map< char, storage::List< storage::Triplet< biol::SSType, ResidueSimple, ResidueSimple> > >::const_iterator chain_itr( SSE_DEFINITIONS.Begin()), chain_itr_end( SSE_DEFINITIONS.End());
        chain_itr != chain_itr_end;
        ++chain_itr
      )
      {
        std::set< storage::Triplet< biol::SSType, ResidueSimple, ResidueSimple>, SSEInfoCompareType> &this_set
        (
          set_of_sse_definitions.insert
          (
            std::make_pair
            (
              chain_itr->first,
              std::set< storage::Triplet< biol::SSType, ResidueSimple, ResidueSimple>, SSEInfoCompareType>( std::ptr_fun( &SSEInformationCompare))
            )
          ).first->second
        );

        for
        (
          storage::List< storage::Triplet< biol::SSType, ResidueSimple, ResidueSimple> >::const_iterator itr( chain_itr->second.Begin()), itr_end( chain_itr->second.End());
          itr != itr_end;
          ++itr
        )
        {
          // create the sse definition
          storage::Triplet< biol::SSType, ResidueSimple, ResidueSimple> sse_definition( *itr);

          // check that residues are given in correct order
          if( sse_definition.Third() < sse_definition.Second())
          {
            BCL_MessageCrt
            (
              "SSE definition might be given in wrong order!\n" + util::Format()( *itr)
            );
          }

          std::pair< std::set< storage::Triplet< biol::SSType, ResidueSimple, ResidueSimple>, SSEInfoCompareType>::iterator, bool>
            insert_itr_success;

          // as long as it is impossible to store the sse definition, there might be an overlapping definition
          // one overlapping residue is permitted by the comparison function, more than one is not
          while( !( insert_itr_success = this_set.insert( sse_definition)).second)
          {
            BCL_MessageCrt
            (
              "similar sse definition was already found in pdb. Tried to insert:\n" + util::Format()( sse_definition) +
              "\n which is similar to already existing:\n" + util::Format()( *insert_itr_success.first) +
              "\n they will be merged into one"
            );

            // check that overlapping definitions are of the same sse type
            if( insert_itr_success.first->First() != sse_definition.First())
            {
              BCL_MessageCrt( "overlapping sse definitions have different sstypes");
              break;
            };

            // if sse definition comes before the one already in the set
            if( SSEInformationCompare( sse_definition, *insert_itr_success.first))
            {
              // the sse definition becomes the combination of both
              sse_definition =
                storage::Triplet< biol::SSType, ResidueSimple, ResidueSimple>
                (
                  sse_definition.First(),
                  sse_definition.Second(),
                  insert_itr_success.first->Third()
                );
            }
            else
            {
              sse_definition =
                storage::Triplet< biol::SSType, ResidueSimple, ResidueSimple>
                (
                  sse_definition.First(),
                  insert_itr_success.first->Second(),
                  sse_definition.Third()
                );
            }

            // remove the overlapping definition already in the set
            this_set.erase( insert_itr_success.first);
          }
        }
      }

      storage::Map< char, storage::List< storage::Triplet< biol::SSType, ResidueSimple, ResidueSimple> > > result;

      // convert temp storage to final result
      for
      (
        std::map< char, std::set< storage::Triplet< biol::SSType, ResidueSimple, ResidueSimple>, SSEInfoCompareType> >::const_iterator
          itr( set_of_sse_definitions.begin()), itr_end( set_of_sse_definitions.end());
        itr != itr_end;
        ++itr
      )
      {
        result[ itr->first] = storage::List< storage::Triplet< biol::SSType, ResidueSimple, ResidueSimple> >( itr->second.begin(), itr->second.end());
      }

      return result;
    }

    //! @brief compare the sse infos for overlap
    //! @param SSE_INFO_LHS left hand side sse info
    //! @param SSE_INFO_RHS right hand side sse info
    bool Head::SSEInformationCompare
    (
      const storage::Triplet< biol::SSType, ResidueSimple, ResidueSimple> &SSE_INFO_LHS,
      const storage::Triplet< biol::SSType, ResidueSimple, ResidueSimple> &SSE_INFO_RHS
    )
    {
      // in the case that two sses share one and the same residue
      if( SSE_INFO_LHS.Third() == SSE_INFO_RHS.Second())
      {
        // compare the initial residue of SSE_INFO_LHS with Terminating residue if SSE_INFO_RHS
        return SSE_INFO_LHS.Second() < SSE_INFO_RHS.Third();
      }

      // return if LHS comes before RHS
      return SSE_INFO_LHS.Third() < SSE_INFO_RHS.Second();
    }

  } // namespace pdb
} // namespace bcl

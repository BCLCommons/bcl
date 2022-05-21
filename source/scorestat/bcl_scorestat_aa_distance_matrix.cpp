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

// include header for this class
#include "scorestat/bcl_scorestat_aa_distance_matrix.h"
// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_ensemble.h"
#include "io/bcl_io_serialization.h"
#include "storage/bcl_storage_table.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace scorestat
  {
    //! single instance of this class
    const util::SiPtr< util::ObjectInterface> AADistanceMatrix::s_Instance
    (
      util::Enumerated< assemble::AnalyzeProteinEnsembleInterface>::AddInstance( new AADistanceMatrix())
    );

    //! @brief OutputOption as string
    //! @param OUTPUT_OPTION the OutputOption
    //! @return the string for the OutputOption
    const std::string &AADistanceMatrix::GetOutputOptionName( const OutputOption &OUTPUT_OPTION)
    {
      static const std::string s_names[] =
      {
          "Table",
          GetStaticClassName< AADistanceMatrix::OutputOption>()
      };
      return s_names[ size_t( OUTPUT_OPTION)];
    }

    //! @brief Output filename as string
    //! @param OutputOption the desired Output Type
    //! @return the string for the output file extension
    const std::string &AADistanceMatrix::GetOutputFileName( const OutputOption &OUTPUT_OPTION)
    {
      static const std::string s_names[] =
      {
        "aa_distance_matrix.tbl",
        GetStaticClassName< AADistanceMatrix::OutputOption>()
      };
      return s_names[ size_t( OUTPUT_OPTION)];
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

      //! @brief implicit constructor
      AADistanceMatrix::AADistanceMatrix() :
          m_OutputOption( e_Table),
          m_ColumnChainId( ""),
          m_RowChainId( "")
      {
        // nothing else to do
      }

      //! @brief virtual copy constructor
      //! @return a pointer to an instance of the copied AADistanceMatrix
      AADistanceMatrix *AADistanceMatrix::Clone() const
      {
        return new AADistanceMatrix( *this);
      }

  /////////////////
  // data access //
  /////////////////

      //! @brief returns class name
      //! @return the class name as constant reference to std::string
      const std::string &AADistanceMatrix::GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief gives the string to append to the the end of a filename to identify this analysis
      //! @return the string to append to the the end of a filename to identify this analysis
      const std::string &AADistanceMatrix::GetOutFilePostfix() const
      {
        return GetOutputFileName( m_OutputOption);
      }

      //! @brief returns chain ids
      //! @return chain ids
      const std::string &AADistanceMatrix::GetColumnChainId() const
      {
        return m_ColumnChainId;
      }

      //! @brief returns chain ids
      //! @return chain ids
      const std::string &AADistanceMatrix::GetRowChainId() const
      {
        return m_RowChainId;
      }

      //! @brief returns the name used for this class in an object data label
      //! @return the name used for this class in an object data label
      const std::string &AADistanceMatrix::GetAlias() const
      {
        static const std::string s_name( "AADistanceMatrix");
        return s_name;
      }

  ///////////////
  // operators //
  ///////////////

      //! @brief takes an ensemble and does an analysis resulting in string containing formatted information
      //! @param ENSEMBLE the protein ensemble that will be analyzed
      //! @return string which has the analyzed information about the ensemble
      std::string AADistanceMatrix::operator()( const assemble::ProteinEnsemble &ENSEMBLE) const
      {
        // make sure the protein ensemble is not empty
        BCL_Assert( !ENSEMBLE.IsEmpty(), "protein ensemble is empty");

        // table header of the matrix for storing aapair distances
        storage::Table< double> distance_matrix;
        storage::Vector< std::string> matrix_column_names;
        std::map< std::pair< std::string, std::string>, double> distance_map;

        // get first protein model from the ensemble
        const assemble::ProteinModel first_protein_model( ( **ENSEMBLE.Begin()));

        // get the desired row chain
        const util::ShPtr< biol::AASequence> row_chain( first_protein_model.GetChain( m_RowChainId[ 0])->GetSequence());

        // get the desired column chain
        const util::ShPtr< biol::AASequence> column_chain( first_protein_model.GetChain( m_ColumnChainId[ 0])->GetSequence());

        // make sure protein model is qualified
        BCL_Assert
        (
          row_chain.IsDefined() && column_chain.IsDefined(),
          "Unqualified protein model: the current protein model does not have desired chains"
        );

        // iterate through the column chain
        for
        (
          util::ShPtrVector< biol::AABase>::const_iterator
            col_aa_itr( column_chain->Begin()), col_aa_itr_end( column_chain->End());
          col_aa_itr != col_aa_itr_end;
          ++col_aa_itr
        )
        {
          // skip undefined residues
          if( !( *col_aa_itr)->GetType()->IsNaturalAminoAcid() || ( *col_aa_itr)->GetAtomCoordinates().IsEmpty())
          {
            continue;
          }
          matrix_column_names.PushBack
          (
            ( *col_aa_itr)->GetType()->GetOneLetterCode() + util::Format()( ( *col_aa_itr)->GetSeqID())
          );
        }

        // initialize the table header
        const util::ShPtr< storage::TableHeader> sp_table_header( new storage::TableHeader( matrix_column_names));
        distance_matrix = storage::Table< double>( sp_table_header);

        // iterate through all protein models in the ensemble
        size_t number_valide_protein_models( 0);
        for
        (
          util::ShPtrVector< assemble::ProteinModel>::const_iterator
            protein_model_itr( ENSEMBLE.Begin()), protein_model_itr_end( ENSEMBLE.End());
          protein_model_itr != protein_model_itr_end;
          ++protein_model_itr
        )
        {
          // get current protein model
          const assemble::ProteinModel protein_model( **protein_model_itr);

          // get the desired row chain
          const util::ShPtr< biol::AASequence> row_chain( protein_model.GetChain( m_RowChainId[ 0])->GetSequence());

          // get the desired column chain
          const util::ShPtr< biol::AASequence> column_chain( protein_model.GetChain( m_ColumnChainId[ 0])->GetSequence());

          // skip unqualified protein models
          if( !row_chain.IsDefined() || !column_chain.IsDefined())
          {
            BCL_MessageStd( "Unqualified protein model: the current protein model does not have desired chains");
            continue;
          }

          // increment the number of valide protein models
          ++number_valide_protein_models;

          // iterate through row residues to insert rows
          for
          (
            util::ShPtrVector< biol::AABase>::const_iterator
              row_aa_itr( row_chain->Begin()), row_aa_itr_end( row_chain->End());
            row_aa_itr != row_aa_itr_end;
            ++row_aa_itr
          )
          {
            // skip undefined residues
            if( !( *row_aa_itr)->GetType()->IsNaturalAminoAcid() || ( *row_aa_itr)->GetAtomCoordinates().IsEmpty())
            {
              continue;
            }

            // row_aa identifier
            const std::string row_aa_id
            (
              ( *row_aa_itr)->GetType()->GetOneLetterCode() + util::Format()( ( *row_aa_itr)->GetSeqID())
            );

            // iterate through column residues
            for
            (
              util::ShPtrVector< biol::AABase>::const_iterator
                col_aa_itr( column_chain->Begin()), col_aa_itr_end( column_chain->End());
              col_aa_itr != col_aa_itr_end;
              ++col_aa_itr
            )
            {
              // skip undefined residues
              if( !( *col_aa_itr)->GetType()->IsNaturalAminoAcid() || ( *col_aa_itr)->GetAtomCoordinates().IsEmpty())
              {
                continue;
              }

              // col_aa identifier
              const std::string col_aa_id
              (
                ( *col_aa_itr)->GetType()->GetOneLetterCode() + util::Format()( ( *col_aa_itr)->GetSeqID())
              );

              distance_map[ std::pair< std::string, std::string>( row_aa_id, col_aa_id)] +=
                  biol::FirstSidechainAtomDistance( ( **row_aa_itr), ( **col_aa_itr));
            }
          }
        }

        // iterate through row residues to insert rows
        for
        (
          util::ShPtrVector< biol::AABase>::const_iterator
            row_aa_itr( row_chain->Begin()), row_aa_itr_end( row_chain->End());
          row_aa_itr != row_aa_itr_end;
          ++row_aa_itr
        )
        {
          // skip undefined residues
          if( !( *row_aa_itr)->GetType()->IsNaturalAminoAcid() || ( *row_aa_itr)->GetAtomCoordinates().IsEmpty())
          {
            continue;
          }

          // row_aa identifier
          const std::string row_aa_id
          (
            ( *row_aa_itr)->GetType()->GetOneLetterCode() + util::Format()( ( *row_aa_itr)->GetSeqID())
          );

          // create vector for holding distances
          storage::Vector< double> distances;

          // iterate through column residues
          for
          (
            util::ShPtrVector< biol::AABase>::const_iterator
              col_aa_itr( column_chain->Begin()), col_aa_itr_end( column_chain->End());
            col_aa_itr != col_aa_itr_end;
            ++col_aa_itr
          )
          {
            // col_aa identifier
            const std::string col_aa_id
            (
              ( *col_aa_itr)->GetType()->GetOneLetterCode() + util::Format()( ( *col_aa_itr)->GetSeqID())
            );

            distances.PushBack
            (
              distance_map[ std::pair< std::string, std::string>( row_aa_id, col_aa_id)] / number_valide_protein_models
            );
          }

          // insert row to table
          distance_matrix.InsertRow
          (
            row_aa_id,
            distances,
            true
          );

        }
        // write distance matrix
        std::ostringstream stream;
        if( m_OutputOption == e_Table)
        {
          distance_matrix.WriteFormatted( stream);
        }

        return stream.str();
      }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer AADistanceMatrix::GetSerializer() const
    {
      io::Serializer paramters;

      paramters.SetClassDescription
      (
        "Compute amino acid distance matrix averaged over a protein conformational ensemble."
      );

      paramters.AddInitializer
      (
        "row_chain_id",
        "id of the chain containing residues corresponding to row names",
        io::Serialization::GetAgent( &m_RowChainId),
        "A"
      );

      paramters.AddInitializer
      (
        "column_chain_id",
        "id of the chain containing residues corresponding to column names",
        io::Serialization::GetAgent( &m_ColumnChainId),
        "A"
      );

      paramters.AddInitializer
      (
        "output",
        "output format",
        io::Serialization::GetAgent( &m_OutputOption),
        "Table"
      );

      return paramters;
    }
  } // namespace scorestat
} // namespace bcl

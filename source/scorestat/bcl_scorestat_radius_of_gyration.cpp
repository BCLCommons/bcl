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
#include "scorestat/bcl_scorestat_radius_of_gyration.h"
// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_ensemble.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_serialization.h"
#include "pdb/bcl_pdb_factory.h"
#include "score/bcl_score_radius_of_gyration.h"
#include "util/bcl_util_wrapper.h"

namespace bcl
{
  namespace scorestat
  {
    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> RadiusOfGyration::s_Instance
    (
      util::Enumerated< assemble::AnalyzeProteinEnsembleInterface>::AddInstance( new RadiusOfGyration())
    );

    //! @brief OutputOption as string
    //! @param OUTPUT_OPTION the OutputOption
    //! @return the string for the OutputOption
    const std::string &RadiusOfGyration::GetOutputOptionName( const OutputOption &OUTPUT_OPTION)
    {
      static std::string s_names[] =
      {
          "Table",
          "Histogram",
          GetStaticClassName< RadiusOfGyration::OutputOption>()
      };
      return s_names[ size_t( OUTPUT_OPTION)];
    }

    //! @brief Output filename as string
    //! @param OutputOption the desired Output Type
    //! @return the string for the output file extension
    const std::string &RadiusOfGyration::GetOutputFileName( const OutputOption &OUTPUT_OPTION)
    {
      static std::string s_output_filename_extenstions[] =
      {
          "radius_of_gyration.tbl",
          "radius_of_gyration.histograms",
          GetStaticClassName< RadiusOfGyration::OutputOption>()
      };
      return s_output_filename_extenstions[ size_t( OUTPUT_OPTION)];
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    RadiusOfGyration::RadiusOfGyration() :
        m_OutputOption( e_Table),
        m_ChainIds( "")
    {
      // nothing to do
    }

    //! @brief virtual copy constructor
    RadiusOfGyration *RadiusOfGyration::Clone() const
    {
      return new RadiusOfGyration( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as constant reference to std::string
    const std::string &RadiusOfGyration::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gives the string to append to the the end of a filename to identify this analysis
    //! @return gives the string to append to the the end of a filename to identify this analysis
    const std::string &RadiusOfGyration::GetOutFilePostfix() const
    {
      return GetOutputFileName( m_OutputOption);
    }

    //! @brief returns chain ids
    //! @return chain ids
    const std::string &RadiusOfGyration::GetChainIds() const
    {
      return m_ChainIds;
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &RadiusOfGyration::GetAlias() const
    {
      static std::string s_Name( "RadiusOfGyration");
      return s_Name;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief takes an ensemble and does an analysis resulting in string containing formatted information
    //! @param ENSEMBLE the protein ensemble that will be analyzed
    //! @return string which has the analyzed information about the ensemble
    std::string RadiusOfGyration::operator()( const assemble::ProteinEnsemble &ENSEMBLE) const
    {
      // make sure the protein ensemble is not empty
      BCL_Assert( !ENSEMBLE.IsEmpty(), "protein ensemble is empty");

      // create table for holding radius of gyration statistics
      storage::Table< double> radius_gyration_table
      (
        storage::TableHeader
        (
          storage::Vector< std::string>::Create
          (
            "subunits", "nr_coordinates", "nr_aa", "rgyr_sqr", "density"
          )
        )
      );

      // iterate over the protein ensemble
      for
      (
        util::ShPtrVector< assemble::ProteinModel>::const_iterator
          protein_model_itr( ENSEMBLE.Begin()), protein_model_itr_end( ENSEMBLE.End());
        protein_model_itr != protein_model_itr_end;
        ++protein_model_itr
      )
      {
        // get the current protein model
        const assemble::ProteinModel &current_protein_model( **protein_model_itr);

        // get current pdb name
        const util::ShPtr< util::Wrapper< std::string> >
          &sp_model_name( current_protein_model.GetProteinModelData()->GetData( assemble::ProteinModelData::e_PDBFile));
        const std::string &model_name( sp_model_name.IsDefined() ? io::File::RemovePath( sp_model_name->GetData()) : "");

        BCL_MessageStd( "Computing Rgyr terms on " + model_name);
        // get all chains of current protein model
        const util::ShPtrVector< assemble::Chain> &all_chains( current_protein_model.GetChains());

        // create a vector of coordinates for defined first side chain atoms of all chains
        util::SiPtrVector< const linal::Vector3D> defined_first_sc_atom_coordinates_all_chains;

        // number of chains used in calculation
        size_t number_of_chains( 0);

        // total number of residues, even outside chains
        size_t number_res( 0);

        // iterate over chains
        for
        (
          util::ShPtrVector< assemble::Chain>::const_iterator chain_itr( all_chains.Begin()), chain_itr_end( all_chains.End());
            chain_itr != chain_itr_end;
          ++chain_itr
        )
        {
          // skip undesired chains, if m_ChainIds us empty, then analyze all chains
          if( !m_ChainIds.empty() && m_ChainIds.find( ( *chain_itr)->GetChainID()) == std::string::npos)
          {
             continue;
          }

          // get coordinate for all side chain atoms
          const util::SiPtrVector< const linal::Vector3D>
            all_first_sc_atom_coordinates_current_chain( ( *chain_itr)->GetAtomCoordinates( biol::GetAtomTypes().GetFirstSidechainAtomTypes()));

          if( all_first_sc_atom_coordinates_current_chain.IsEmpty())
          {
            continue;
          }

          auto sses( ( *chain_itr)->GetSSEs());
          number_res += sses.LastElement()->GetFirstAA()->GetSeqID() - sses.FirstElement()->GetFirstAA()->GetSeqID() + 1;

          // create a vector of coordinates for defined first side chain atoms on current chain
          util::SiPtrVector< const linal::Vector3D> defined_first_sc_atom_coordinates_current_chain;

          // push back the defined ones into the vector
          for
          (
            util::SiPtrVector< const linal::Vector3D>::const_iterator
              coordinate_itr( all_first_sc_atom_coordinates_current_chain.Begin()),
              coordinate_itr_end( all_first_sc_atom_coordinates_current_chain.End());
            coordinate_itr != coordinate_itr_end;
            ++coordinate_itr
          )
          {
            if( ( *coordinate_itr)->IsDefined())
            {
              defined_first_sc_atom_coordinates_current_chain.PushBack( *coordinate_itr);
            }
          }

          // add the vector for the current chain to the vector for all chains
          defined_first_sc_atom_coordinates_all_chains.Append( defined_first_sc_atom_coordinates_current_chain);

          // increment number of chains used in calculation
          ++number_of_chains;
        }

        // create variable for holding radius of gyration of current model
        double radius_of_gyration( 0.0);

        // if membrane flag is given
        if( biol::Membrane::GetFlagMembrane()->GetFlag())
        {
          // calculate the collapsed radius of gyration
          radius_of_gyration = score::RadiusOfGyration::SquareRadiusOfGyrationCollapsed
              (
                defined_first_sc_atom_coordinates_all_chains,
                current_protein_model.GetProteinModelData()->GetData( assemble::ProteinModelData::e_Membrane)
              );
        }
        // else soluble proteins
        else
        {
          // calculate radius of gyration in the standard way
          radius_of_gyration = coord::SquareRadiusOfGyration( defined_first_sc_atom_coordinates_all_chains);
        }

        const double atom_vdw_radius( 1.5);

        const double est_volume
        (
          std::max
          (
            std::max
            (
              coord::EstimateVolume( defined_first_sc_atom_coordinates_all_chains, 3.4, 10.0, atom_vdw_radius),
              1.0
            ),
            4.0 / 3.0 * math::g_Pi / 0.74 * number_res * atom_vdw_radius * atom_vdw_radius * atom_vdw_radius
          )
        );

        // insert into table
        radius_gyration_table.InsertRow
        (
          model_name,
          storage::Vector< double>::Create
          (
            number_of_chains,
            defined_first_sc_atom_coordinates_all_chains.GetSize(),
            number_res,
            radius_of_gyration,
            double( number_res) / est_volume
          )
        );
      }

      // write statistics
      std::ostringstream stream;

      // if output option is table
      if( m_OutputOption == e_Table)
      {
        radius_gyration_table.WriteFormatted( stream);
      }
      else if( m_OutputOption == e_Histogram)
      {
        // create histograms that holds radius of gyration statistics
        const storage::Map< std::string, math::Histogram> radius_of_gyration_histograms
        (
          score::RadiusOfGyration::HistogramsFromTable( radius_gyration_table)
        );

        // iterate over histograms
        for
        (
          storage::Map< std::string, math::Histogram>::const_iterator
            histogram_itr( radius_of_gyration_histograms.Begin()), histogram_itr_end( radius_of_gyration_histograms.End());
          histogram_itr != histogram_itr_end;
          ++histogram_itr
        )
        {
          stream << histogram_itr->first + " histogram" << '\n';
          stream << histogram_itr->second << '\n';
        }
      }

      return stream.str();
    } // end of operator ()

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer RadiusOfGyration::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Computes radius of gyration statistics."
      );

      parameters.AddInitializer
      (
        "chain_ids",
        "string of chain ids to be analyzed, if not specified, take all chains",
        io::Serialization::GetAgent( &m_ChainIds),
        ""
      );

      parameters.AddInitializer
      (
        "output",
        "what type of outputs to provide",
        io::Serialization::GetAgent( &m_OutputOption),
        "Table"
      );
      return parameters;
    }
  } // namespace scorestat
} // namespace bcl


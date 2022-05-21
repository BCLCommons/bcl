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
#include "sspred/bcl_sspred_pdb.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "biol/bcl_biol_atom.h"
#include "biol/bcl_biol_membrane.h"
#include "command/bcl_command_parameter_interface.h"
#include "io/bcl_io_directory_entry.h"
#include "io/bcl_io_file.h"
#include "util/bcl_util_wrapper.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace sspred
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> PDB::s_Instance
    (
      GetObjectInstances().AddInstance( new PDB())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    PDB::PDB() :
      m_SSType( biol::GetSSTypes().e_Undefined),
      m_EnvironmentType( biol::GetEnvironmentTypes().e_Undefined)
    {
    }

    //! @brief constructor from a given SSType enum and an EnvironmentType enum
    //! @param SS_TYPE SSType enum of interest
    //! @param ENVIRONMENT_TYPE SSType EnvironmentType enum
    PDB::PDB( const biol::SSType &SS_TYPE, const biol::EnvironmentType &ENVIRONMENT_TYPE) :
      m_SSType( SS_TYPE),
      m_EnvironmentType( ENVIRONMENT_TYPE)
    {
    }

    //! @brief Clone function
    //! @return pointer to new PDB
    PDB *PDB::Clone() const
    {
      return new PDB( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &PDB::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief get three state, environment independent secondary structure prediction
    //! @return three state, environment independent secondary structure prediction
    linal::Vector3D PDB::GetThreeStatePrediction() const
    {
      return m_SSType->GetThreeStatePrediction();
    }

    //! @brief get nine state secondary structure prediction ( 3 SSTypes for all 3 EnvironmentTypes)
    //! @return three state secondary structure prediction
    linal::Matrix< double> PDB::GetNineStatePrediction() const
    {
      return ConvertThreeStateToNineState( GetThreeStatePrediction(), m_EnvironmentType);
    }

    //! @brief find the SSType with highest prediction and returns it
    //! @return SSType with highest prediction
    biol::SSType PDB::GetOneStateSSPrediction() const
    {
      return m_SSType;
    }

    //! @brief find the TMTypes with highest prediction and returns it
    //! @return TMTYpe with highest prediction
    biol::EnvironmentType PDB::GetOneStateTMPrediction() const
    {
      return m_EnvironmentType;
    }

    //! @brief find the SSType-TMType pair with highest prediction and returns it
    //! @return SSType-TMType pair with highest prediction
    storage::Pair< biol::SSType, biol::EnvironmentType> PDB::GetOneStateSSTMPrediction() const
    {
      return storage::Pair< biol::SSType, biol::EnvironmentType>( m_SSType, m_EnvironmentType);
    }

    //! @brief set the environment types for the residues in the protein model
    //! @param PROTEIN_MODEL model containing residues to be set
    //! @param USE_PDBTM_MEMBRANE_THICKNESS true to use membrane thickness from the PDBTM if it was read in
    //! @return true if the membrane environment was found; false otherwise
    bool PDB::SetEnvironmentTypes( assemble::ProteinModel &PROTEIN_MODEL, const bool &USE_PDBTM_MEMBRANE_THICKNESS)
    {
      // get membrane
      util::ShPtr< biol::Membrane> sp_membrane
      (
        PROTEIN_MODEL.GetProteinModelData()->GetData( assemble::ProteinModelData::e_Membrane)
      );

      // if no membrane was found, try to load it
      if( !sp_membrane.IsDefined())
      {
        // get the filename
        util::ShPtr< util::Wrapper< std::string> > sp_filename_wrapper
        (
          PROTEIN_MODEL.GetProteinModelData()->GetData( assemble::ProteinModelData::e_PDBFile)
        );

        // if the filename was available
        if( sp_filename_wrapper.IsDefined())
        {
          std::string pdb_filename( sp_filename_wrapper->GetData());

          // Remove the last extension
          std::string basename( io::File::RemoveLastExtension( io::File::RemoveCompressionExtension( pdb_filename)));

          // determine the xml name
          storage::Vector< std::string> candidate_xml_names;
          candidate_xml_names.PushBack( basename + ".xml");

          // if the pdb contains only a single chain, try removing the chain id from the file and looking for it
          if( PROTEIN_MODEL.GetChains().GetSize() <= size_t( 1))
          {
            candidate_xml_names.PushBack( basename.substr( 0, basename.size() - 1) + ".xml");
          }

          // also only keeping the first 4 characters of the id, in case the user provided something like 1ubibiobcl.pdb
          storage::VectorND< 2, std::string> pdb_path_basename( io::File::SplitToPathAndFileName( basename));
          if( pdb_path_basename( 1).size() > size_t( 7))
          {
            candidate_xml_names.PushBack
            (
              ( pdb_path_basename.First().empty() ? std::string() : pdb_path_basename.First() + PATH_SEPARATOR) + pdb_path_basename.Second().substr( 0, 4) + ".xml"
            );
          }

          storage::Vector< std::string>::const_iterator itr_filename( candidate_xml_names.Begin());
          // find the first existing xml file
          for
          (
            storage::Vector< std::string>::const_iterator itr_filename_end( candidate_xml_names.End());
            itr_filename != itr_filename_end;
            ++itr_filename
          )
          {
            if( io::DirectoryEntry( *itr_filename).DoesExist())
            {
              break;
            }
          }

          if( itr_filename == candidate_xml_names.End())
          {
            // xml file not found, no membrane to set
            BCL_MessageVrb
            (
              "No membrane XML file found; looked at: " + util::Format()( candidate_xml_names)
            );
            return false;
          }

          // open the xml file
          io::IFStream input;
          io::File::MustOpenIFStream( input, *itr_filename);
          const storage::Pair< biol::Membrane, math::TransformationMatrix3D>
            membrane_transformation
            (
              biol::Membrane::MembraneAndTransformationFromPDBTMXML
              (
                input,
                biol::Membrane::GetParameterTransitionThickness()->GetNumericalValue< double>(),
                biol::Membrane::GetParameterGapThickness()->GetNumericalValue< double>()
              )
            );

          io::File::CloseClearFStream( input);
          if( !membrane_transformation.First().IsDefined())
          {
            // if no membrane was found, just return
            return false;
          }

          PROTEIN_MODEL.ConnectSSEToChainData();

          // transform model with the membrane transformation matrix
          PROTEIN_MODEL.Transform( membrane_transformation.Second());

          // transform all internal sequences

          // loop over all chains and transform the internal sequences as well
          for
          (
            util::ShPtrVector< assemble::Chain>::const_iterator
              chain_itr( PROTEIN_MODEL.GetChains().Begin()), chain_itr_end( PROTEIN_MODEL.GetChains().End());
            chain_itr != chain_itr_end;
            ++chain_itr
          )
          {
            util::ShPtr< assemble::Chain> chain( *chain_itr);
            chain->GetSequence()->Transform( membrane_transformation.Second());
          }

          // override membrane thickness with the constant thickness given on the command line, if it was set
          if( !USE_PDBTM_MEMBRANE_THICKNESS)
          {
            const biol::Membrane &membrane_pdbtm( membrane_transformation.First());
            sp_membrane =
              util::ShPtr< biol::Membrane>
              (
                new biol::Membrane
                (
                  biol::Membrane::GetCommandLineMembrane().GetThicknesses(),
                  membrane_pdbtm.GetNormal()
                )
              );
          }
          else
          {
            // copy the membrane
            sp_membrane = util::ShPtr< biol::Membrane>( membrane_transformation.First().Clone());
          }

          // set membrane
          util::ShPtr< assemble::ProteinModelData> sp_data( PROTEIN_MODEL.GetProteinModelData());
          sp_data->Insert( assemble::ProteinModelData::e_Membrane, sp_membrane);
          PROTEIN_MODEL.SetProteinModelData( sp_data);
        }
      }

      // end if no membrane found
      if( !sp_membrane.IsDefined() || !sp_membrane->IsDefined())
      {
        return false;
      }

      // iterate through the chains
      for
      (
        util::ShPtrVector< assemble::Chain>::const_iterator chain_itr( PROTEIN_MODEL.GetChains().Begin()),
          chain_itr_end( PROTEIN_MODEL.GetChains().End());
        chain_itr != chain_itr_end; ++chain_itr
      )
      {
        // iterate through the sequence
        for
        (
          biol::AASequence::const_iterator aa_itr( ( *chain_itr)->GetSequence()->Begin()),
            aa_itr_end( ( *chain_itr)->GetSequence()->End());
          aa_itr != aa_itr_end; ++aa_itr
        )
        {
          // create a pointer to the data
          util::ShPtr< biol::AABase> sp_aa( *aa_itr);
          const util::SiPtr< const MethodInterface> sp_aa_ss_pdb( sp_aa->GetSSPrediction( GetMethods().e_PDB));
          biol::SSType existing_ss_type( biol::GetSSTypes().COIL);
          if( !sp_aa_ss_pdb.IsDefined())
          {
            BCL_MessageCrt
            (
              "Warning: AA with id: " + sp_aa->GetIdentification() + " had no secondary structure from the pdb file, "
              "setting it to type coil"
            );
          }
          else
          {
            existing_ss_type = sp_aa_ss_pdb->GetOneStateSSPrediction();
          }

          // determine the environment type, given the coordinate
          const linal::Vector3D atom_coord( sp_aa->GetCA().GetCoordinates());

          biol::EnvironmentType env_type( sp_membrane->DetermineEnvironmentType( atom_coord));

          if( !env_type.IsDefined())
          {
            if( atom_coord.IsDefined())
            {
              BCL_MessageCrt
              (
                "Warning: AA with id: " + sp_aa->GetIdentification()
                + " has an undefined membrane environment; z-coordinate was: " + util::Format()( atom_coord)
              );
            }
            env_type = biol::GetEnvironmentTypes().e_Solution;
          }

          // set the ssprediction
          sp_aa->SetSSPrediction( GetMethods().e_PDB, PDB( existing_ss_type, env_type));
        }
      }
      return true;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read secondary structure predictions for given amino acid from the provided ISTREAM
    //! @param ISTREAM input stream
    //! @param AMINO_ACID amino acid into which sspredictions will be read
    //! @return std::istream which was read from
    std::istream &PDB::ReadPredictionsForAA
    (
      std::istream &ISTREAM,
      biol::AABase &AMINO_ACID
    ) const
    {
      return ISTREAM;
    }

    //! @brief read secondary structure predictions for given sequence from the provided ISTREAM
    //! @param ISTREAM input stream
    //! @param AA_SEQUENCE AASequence into which sspredictions will be read
    //! @return std::istream which was read from
    std::istream &PDB::ReadPredictionsForAASequence
    (
      std::istream &ISTREAM,
      biol::AASequence &AA_SEQUENCE
    ) const
    {
      return ISTREAM;
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &PDB::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_SSType, ISTREAM);
      io::Serialize::Read( m_EnvironmentType, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &PDB::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_SSType, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_EnvironmentType, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  } // namespace sspred
} // namespace bcl

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
#include "restraint/bcl_restraint_analyze_per_residue_rmsd.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_locator_atom.h"
#include "assemble/bcl_assemble_protein_ensemble.h"
#include "assemble/bcl_assemble_quality.h"
#include "command/bcl_command_parameter_check_file_existence.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_serialization.h"
#include "io/bcl_io_serialization_via_static_functions.h"
#include "pdb/bcl_pdb_factory.h"
#include "quality/bcl_quality_rmsd.h"
#include "util/bcl_util_wrapper.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> AnalyzePerResidueRMSD::s_Instance
    (
      util::Enumerated< assemble::AnalyzeProteinEnsembleInterface>::AddInstance( new AnalyzePerResidueRMSD())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    AnalyzePerResidueRMSD::AnalyzePerResidueRMSD() :
      m_OutFilePostFix( ".AnalyzePerResidueRMSD"),
      m_SuperimposeMeasure(),
      m_TemplateModel(),
      m_QualityMeasure( quality::GetMeasures().e_RMSD),
      m_PymolOuputFilename( "AnalyzePerResidueRMSD.py")
    {
    }

    //! @brief Clone function
    //! @return pointer to new AnalyzePerResidueRMSD
    AnalyzePerResidueRMSD *AnalyzePerResidueRMSD::Clone() const
    {
      return new AnalyzePerResidueRMSD( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AnalyzePerResidueRMSD::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gives the string to append to the the end of a filename to identify this analysis
    //! @return gives the string to append to the the end of a filename to identify this analysis
    const std::string &AnalyzePerResidueRMSD::GetOutFilePostfix() const
    {
      return m_OutFilePostFix;
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &AnalyzePerResidueRMSD::GetAlias() const
    {
      static const std::string s_name( "PerResidueRMSD");
      return s_name;
    }

  ////////////////
  // operations //
  ////////////////

  ///////////////
  // operators //
  ///////////////

    //! @brief takes an ensemble and does an analysis resulting in string containing formatted information
    //! @param ENSEMBLE the ensemble that will be analyzed
    //! @return string which has the analyzed information about the ensemble
    std::string AnalyzePerResidueRMSD::operator()( const assemble::ProteinEnsemble &ENSEMBLE) const
    {
      // string to hold analysis
      std::string analysis;

      // open write to pymol output script
      io::OFStream write;
      io::File::MustOpenOFStream
      (
        write, GetFlagOutFilePrefix()->GetFirstParameter()->GetValue() + m_PymolOuputFilename
      );

      // write python script loading of template model
      write << "cmd.load( \""
            <<  util::ShPtr< util::Wrapper< std::string> >
                (
                  m_TemplateModel.GetProteinModelData()->GetData( assemble::ProteinModelData::e_PDBFile)
                )->GetData()
            << "\", \"template_model\")\n";

      // map to hold for each residue the corresponding coordinates
      storage::Map< assemble::LocatorAtom, storage::Vector< linal::Vector3D> > atom_coordinate_list
      (
        GetAtomCoordinates( ENSEMBLE, m_SuperimposeMeasure, m_TemplateModel)
      );

      // iterate over the map of coordinates for each residue
      for
      (
        storage::Map< assemble::LocatorAtom, storage::Vector< linal::Vector3D> >::const_iterator
          map_itr( atom_coordinate_list.Begin()), map_itr_end( atom_coordinate_list.End());
        map_itr != map_itr_end;
        ++map_itr
      )
      {
        // current atom
        const assemble::LocatorAtom &atom( map_itr->first);

        // vector of coordinates for the current atom over the whole ensemble
        storage::Vector< linal::Vector3D> coords( map_itr->second);

        BCL_MessageDbg( " atom " + atom.GetIdentification());

        // calculate the rmsd of the coordinates to one another
        const util::SiPtrVector< const linal::Vector3D> all_coords( util::ConvertToSiPtrVector( coords));
        const storage::Pair< storage::VectorND< 2, double>, math::RunningAverageSD< double> > rmsd_info
        (
          quality::RMSD::RealSpaceRMSDPairwise( all_coords)
        );

        // if rmsd is not defined
        if( !util::IsDefined( rmsd_info.First().First()))
        {
          // skip to next
          continue;
        }

        // write the analysis information
        analysis += atom.GetIdentification() + " rmsd " + util::Format().W( 9).FFP( 3)( rmsd_info.First().First()) +
                    " rmsd_stddev " + util::Format().W( 9).FFP( 3)( rmsd_info.First().Second()) +
                    " mean_dist " + util::Format().W( 9).FFP( 3)( rmsd_info.Second().GetAverage()) +
                    " dist_stddev " + util::Format().W( 9).FFP( 3)( rmsd_info.Second().GetStandardDeviation()) + '\n';

        // write to the python script for visualization
        write << "cmd.alter( \"template_model and " << atom.GetPymolResidueSelection() << "\", \"b= "
              << util::Format()( rmsd_info.First().First()) << "\")\n";

      }

      // write the ending information to the python script for visualization
      write <<  "cmd.show_as(\"cartoon\"   ,\"template_model\")\n";
      write <<  "cmd.cartoon(\"putty\"   ,\"template_model\")\n";
      write <<  "cmd.unset(\"cartoon_smooth_loops\"   ,\"template_model\")\n";
      write <<  "cmd.unset(\"cartoon_flat_sheets\"   ,\"template_model\")\n";

      return analysis;
    }

  //////////////////////
  // input and output //
  //////////////////////

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer AnalyzePerResidueRMSD::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "For an ensemble of models, superimposes them onto a template structure using the provided superimpose method. "
        "The rmsd of each residue between the ensemble of models is then calculated. Also outputs a python script that"
        " can be used with pymol to visualize the rmsd over the structure by using the putty representation. Residues"
        " with larger rmsd will be thicker. The template model is used in the pymol script."
      );

      parameters.AddInitializer
      (
        "filename_postfix",
        "the postfix that will be appended to the output file in order to identify this analysis",
        io::Serialization::GetAgent( &m_OutFilePostFix),
        ".AnalyzeAtomDistanceScore"
      );

      parameters.AddInitializer
      (
        "superimpose_method",
        "the method that should be used for superimposing onto the template",
        io::Serialization::GetAgent( &m_SuperimposeMeasure),
        quality::GetSuperimposeMeasures().e_RMSD.GetName()
      );

      parameters.AddInitializer
      (
        "template_pdb_filename",
        "the pdb filename of the protein model the other models will be superimposed onto",
        util::OwnPtr< io::SerializationInterface>
        (
          new io::SerializationViaStaticFunctions< assemble::ProteinModel>
          (
            command::ParameterCheckFileExistence(),
            ( &ProteinModelAsString),
            ( &ProteinModelFromString),
            &m_TemplateModel
          )
        ),
        "template.pdb"
      );

      parameters.AddInitializer
      (
        "pymol_output_filename",
        "The filename of the pymol script that will be outputted showing the accessibility changes",
        io::Serialization::GetAgent( &m_PymolOuputFilename),
        "accessibilities.pml"
      );

      return parameters;
    }

    //! @brief converts a protein model into a string by returning the pdb filename stored in the protein model data
    //! @param MODEL the model that will be converted into a string
    //! @return the string of the pdb filename stored in the protein model data
    std::string AnalyzePerResidueRMSD::ProteinModelAsString( const assemble::ProteinModel &MODEL)
    {
      util::ShPtr< util::Wrapper< std::string> > pdb
      (
        MODEL.GetProteinModelData()->GetData( assemble::ProteinModelData::e_PDBFile)
      );

      if( !pdb.IsDefined())
      {
        return "none.pdb";
      }

      return *pdb;
    }

    //! @brief converts a string into a protein model by assuming the string is a pdb filename
    //! @param MODEL protein model to setup
    //! @param NAME string which is the pdb filename
    //! @param ERR_STREAM stream to write out erros to
    bool AnalyzePerResidueRMSD::ProteinModelFromString
    (
      assemble::ProteinModel &MODEL, const std::string &NAME, std::ostream &ERR_STREAM
    )
    {
      // instantiate the pdb factory
      pdb::Factory factory;

      // instantiate protein model
      MODEL = factory.ProteinModelFromPDBFilename( NAME);

      // set the pdb filename in the protein model data
      util::ShPtr< assemble::ProteinModelData>( MODEL.GetProteinModelData())->Insert
      (
        assemble::ProteinModelData::e_PDBFile,
        util::ShPtr< util::Wrapper< std::string> >( new util::Wrapper< std::string>( NAME))
      );
      return true;
    }

    //! @brief gives the CA atom coordinates for each residue for each model in an ensemble
    //!        Each model in the ensemble is superimposed onto a templete model according to the provided
    //!        superimposition method.
    //! @param ENSEMBLE the models whose ca coordinate will be gotten for each residue
    //! @param SUPERIMPOSE_MEASURE the method of superimposition used for the models of the ensemble onto the template
    //! @param TEMPLATE_MODEL the model that the ensemble models will be superimposed on
    //! @return map with a locator for each atom and the vector of associated coordinates coming from the models in
    //!         the ensemble
    storage::Map< assemble::LocatorAtom, storage::Vector< linal::Vector3D> > AnalyzePerResidueRMSD::GetAtomCoordinates
    (
      const assemble::ProteinEnsemble &ENSEMBLE,
      const quality::SuperimposeMeasure &SUPERIMPOSE_MEASURE,
      const assemble::ProteinModel &TEMPLATE_MODEL
    )
    {

      // map to hold for each residue the corresponding coordinates
      storage::Map< assemble::LocatorAtom, storage::Vector< linal::Vector3D> > atom_coordinate_list;

      // iterate through the ensemble
      for
      (
        assemble::ProteinEnsemble::const_iterator ensemble_itr( ENSEMBLE.Begin()), ensemble_itr_end( ENSEMBLE.End());
        ensemble_itr != ensemble_itr_end;
        ++ensemble_itr
      )
      {
        // copy the protein model
        assemble::ProteinModel model( **ensemble_itr);

        // atom type for superimposing
        static const biol::AtomType s_atom_type( biol::GetAtomTypes().CA);

        // atom types for superimposing
        static const storage::Set< biol::AtomType> s_atom_types( s_atom_type);

        // superimpose the current model onto the template
        assemble::Quality::SuperimposeModel( SUPERIMPOSE_MEASURE, model, TEMPLATE_MODEL, s_atom_types);

        // get the residues from the superimposed model
        const util::SiPtrVector< const biol::AABase> aas( model.GetAminoAcids());

        // iterate through the residues of the current model
        for
        (
          util::SiPtrVector< const biol::AABase>::const_iterator aa_itr( aas.Begin()), aa_itr_end( aas.End());
          aa_itr != aa_itr_end;
          ++aa_itr
        )
        {
          // current residue
          const biol::AABase &current_aa( **aa_itr);

          // create locator atom for current residue
          const assemble::LocatorAtom locator( current_aa.GetChainID(), current_aa.GetSeqID(), s_atom_type);

          // get the coordinates for the atom types
          const util::SiPtrVector< const linal::Vector3D> coords( current_aa.GetAtomCoordinates( s_atom_types));

          // assert there is only one set of coordinates, otherwise calculating rmsd is worse
          BCL_Assert
          (
            coords.GetSize() == 1, "number of coordinates is " + util::Format()( coords.GetSize()) +
            " but probably should be 1"
          );

          if( !coords.FirstElement()->IsDefined())
          {
            continue;
          }

          // add coordinates for the atom into the map
          atom_coordinate_list[ locator].Append
          (
            util::ConvertToStorageVector< linal::Vector3D>( current_aa.GetAtomCoordinates( s_atom_types))
          );
        }
      }

      return atom_coordinate_list;
    }

  } // namespace restraint
} // namespace bcl

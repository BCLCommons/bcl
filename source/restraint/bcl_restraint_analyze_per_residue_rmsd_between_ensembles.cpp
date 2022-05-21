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
#include "restraint/bcl_restraint_analyze_per_residue_rmsd_between_ensembles.h"

// includes from bcl - sorted alphabetically
#include "align/bcl_align_alignment_interface.h"
#include "assemble/bcl_assemble_locator_atom.h"
#include "assemble/bcl_assemble_quality.h"
#include "command/bcl_command_flag_interface.h"
#include "command/bcl_command_parameter_check_file_existence.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_serialization.h"
#include "io/bcl_io_serialization_via_static_functions.h"
#include "quality/bcl_quality_rmsd.h"
#include "restraint/bcl_restraint_analyze_accessibility_change.h"
#include "restraint/bcl_restraint_analyze_per_residue_rmsd.h"
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
    const util::SiPtr< const util::ObjectInterface> AnalyzePerResidueRMSDBetweenEnsembles::s_Instance
    (
      util::Enumerated< assemble::AnalyzeProteinEnsembleInterface>::AddInstance( new AnalyzePerResidueRMSDBetweenEnsembles())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    AnalyzePerResidueRMSDBetweenEnsembles::AnalyzePerResidueRMSDBetweenEnsembles() :
      m_OutFilePostFix( ".AnalyzePerResidueRMSDBetweenEnsembles"),
      m_SuperimposeMeasure(),
      m_TemplateModel(),
      m_QualityMeasure( quality::GetMeasures().e_RMSD),
      m_PymolOuputFilename( "AnalyzePerResidueRMSDBetweenEnsembles.py"),
      m_StartEnsemble(),
      m_NormalizeByInternalRMSD( false)
    {
    }

    //! @brief Clone function
    //! @return pointer to new AnalyzePerResidueRMSDBetweenEnsembles
    AnalyzePerResidueRMSDBetweenEnsembles *AnalyzePerResidueRMSDBetweenEnsembles::Clone() const
    {
      return new AnalyzePerResidueRMSDBetweenEnsembles( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AnalyzePerResidueRMSDBetweenEnsembles::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gives the string to append to the the end of a filename to identify this analysis
    //! @return gives the string to append to the the end of a filename to identify this analysis
    const std::string &AnalyzePerResidueRMSDBetweenEnsembles::GetOutFilePostfix() const
    {
      return m_OutFilePostFix;
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &AnalyzePerResidueRMSDBetweenEnsembles::GetAlias() const
    {
      static const std::string s_name( "PerResidueRMSDBetweenEnsembles");
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
    std::string AnalyzePerResidueRMSDBetweenEnsembles::operator()( const assemble::ProteinEnsemble &ENSEMBLE) const
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
      storage::Map< assemble::LocatorAtom, math::RunningAverageSD< double> > atom_coordinate_list;

      // for holding atom lists if want to normalize by the internal per residue rmsd of ENSEMBLE
      storage::Map< assemble::LocatorAtom, storage::Vector< linal::Vector3D> > atom_coordinate_list_norm;
      if( m_NormalizeByInternalRMSD)
      {
        // map to hold for each residue the corresponding coordinates
        atom_coordinate_list_norm = AnalyzePerResidueRMSD::GetAtomCoordinates( ENSEMBLE, m_SuperimposeMeasure, m_TemplateModel);
      }

      // iterate through the current ensemble
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

        // iterate over the start ensemble
        for
        (
          assemble::ProteinEnsemble::const_iterator
            start_ensemble_itr( m_StartEnsemble.Begin()), start_ensemble_itr_end( m_StartEnsemble.End());
            start_ensemble_itr != start_ensemble_itr_end;
          ++start_ensemble_itr
        )
        {
          // superimpose the current model onto the start structure
          assemble::Quality::SuperimposeModel( m_SuperimposeMeasure, model, **start_ensemble_itr, s_atom_types);

          // get alignments for each chain
          storage::Map< char, util::ShPtr< align::AlignmentInterface< biol::AABase> > > alignments
          (
            assemble::Quality::CreateAlignmentProteinModels( model, **start_ensemble_itr)
          );

          // iterate over the chain alignments
          for
          (
            storage::Map< char, util::ShPtr< align::AlignmentInterface< biol::AABase> > >::const_iterator
              chain_itr( alignments.Begin()), chain_itr_end( alignments.End());
            chain_itr != chain_itr_end;
            ++chain_itr
          )
          {
            util::ShPtrList< align::Assignment< biol::AABase> > assignments( chain_itr->second->GetAssignments());

            // iterate over the assignments
            for
            (
              util::ShPtrList< align::Assignment< biol::AABase> >::const_iterator
                assignment_itr( assignments.Begin()), assignment_itr_end( assignments.End());
              assignment_itr != assignment_itr_end;
              ++assignment_itr
            )
            {
              util::SiPtrList< const biol::AABase> members( ( *assignment_itr)->GetMembers());

              // to hold the coords in this assignment
              util::SiPtrVector< const linal::Vector3D> coords;

              // create locator atom for current residue
              assemble::LocatorAtom locator;

              // iterate over the members to get the coordinates
              for
              (
                util::SiPtrList< const biol::AABase>::const_iterator
                  member_itr( members.Begin()), member_itr_end( members.End());
                member_itr != member_itr_end;
                ++member_itr
              )
              {
                const util::SiPtr< const biol::AABase> &aa_base( *member_itr);

                // if the residue is defined
                if( aa_base.IsDefined())
                {
                  // add coordinates to list of coordinates
                  util::SiPtrVector< const linal::Vector3D> current_coord( aa_base->GetAtomCoordinates( s_atom_types));

                  // make sure one coordinate was gotten
                  BCL_Assert( current_coord.GetSize() == 1, "cord size not 1");
                  const linal::Vector3D &ca_coord( *current_coord.FirstElement());

                  // true if the coordinates of the atom are defined
                  if( ca_coord.IsDefined())
                  {
                    // add coords and set locator (locator will be set with each iteration but will be the same)
                    coords.PushBack( ca_coord);
                    locator = assemble::LocatorAtom( aa_base->GetChainID(), aa_base->GetSeqID(), s_atom_type);
                  }
                }
              }

              // should be two coordinates in the list, one for atom from start model one from end model
              if( coords.GetSize() == 2)
              {
                // get the rmsd
                const storage::Pair< storage::VectorND< 2, double>, math::RunningAverageSD< double> > rmsd_info
                (
                  quality::RMSD::RealSpaceRMSDPairwise( coords)
                );

                // for subtracting the background that is the variation within the ENSEMBLE for this residue
                double normalization_factor( 0);
                if( m_NormalizeByInternalRMSD)
                {
                  // try to find the current residue
                  storage::Map< assemble::LocatorAtom, storage::Vector< linal::Vector3D> >::const_iterator norm_itr
                  (
                    atom_coordinate_list_norm.Find( locator)
                  );

                  // true if the residue was found
                  if( norm_itr != atom_coordinate_list_norm.End())
                  {
                    // calculate the rmsd of the coordinates to one another within ENSEMBLE
                    storage::Vector< linal::Vector3D> coords( norm_itr->second);
                    const util::SiPtrVector< const linal::Vector3D> all_coords( util::ConvertToSiPtrVector( coords));
                    const storage::Pair< storage::VectorND< 2, double>, math::RunningAverageSD< double> > rmsd_info_norm
                    (
                      quality::RMSD::RealSpaceRMSDPairwise( all_coords)
                    );
                    // rmsd is the background normalization factor
                    normalization_factor = rmsd_info_norm.First().First();
                  }
                }

                // current rmsd between the ensembles is the actual rmsd between the two ensembles minus the
                // normalization factor
                const double current_rmsd( rmsd_info.First().First() - normalization_factor);

                BCL_MessageDbg
                (
                  "considering rmsd " + util::Format()( current_rmsd) +
                   " for residue " + locator.GetIdentification()
                );

                // add the rmsd to the map
                atom_coordinate_list[ locator] += current_rmsd;
              }
              else //< could be that coordinate is not defined for either of the models, then can't calculate rmsd
              {
                continue;
              }
            } // assignments
          } // alignment chain map
        } // ensemble
      } // ensemble

      // iterate over the map of rmsds for each residue
      for
      (
        storage::Map< assemble::LocatorAtom, math::RunningAverageSD< double> >::const_iterator
          map_itr( atom_coordinate_list.Begin()), map_itr_end( atom_coordinate_list.End());
        map_itr != map_itr_end;
        ++map_itr
      )
      {
        // current atom
        const assemble::LocatorAtom &atom( map_itr->first);

        // write the analysis information
        analysis += atom.GetIdentification() + " rmsd " + util::Format().W( 9).FFP( 3)( map_itr->second.GetAverage()) +
                    " rmsd_stddev " + util::Format().W( 9).FFP( 3)( map_itr->second.GetStandardDeviation()) + '\n';

        // write to the python script for visualization
        write << "cmd.alter( \"template_model and " << atom.GetPymolResidueSelection() << "\", \"b= "
              << util::Format()( map_itr->second.GetAverage()) << "\")\n";
      }

      // write the ending information to the python script for visualization
      write <<  "cmd.show_as(\"cartoon\"   ,\"template_model\")\n";
      write <<  "cmd.cartoon(\"putty\"   ,\"template_model\")\n";
      write <<  "cmd.unset(\"cartoon_smooth_loops\"   ,\"template_model\")\n";
      write <<  "cmd.unset(\"cartoon_flat_sheets\"   ,\"template_model\")\n";

      // return the analysis string
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
    io::Serializer AnalyzePerResidueRMSDBetweenEnsembles::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Calculates the per residue RMSD between two sets of ensembles. CA atoms are used in the calculation. Each "
        " model in each ensemble is pairwise superimposed according to the provided superimposition method, then the "
        " per residue rmsd is calculated along the sequence. These rmsds are averaged over all pairwise RMSDs between "
        " the ensembles. Also outputs a python file that can be opened with pymol in order to vizualize the regions "
        " of the provided template structure that have large per residue RMSDs. This is depicted as thickness using "
        " the putty representation in pymol."
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
            ( &AnalyzePerResidueRMSD::ProteinModelAsString),
            ( &AnalyzePerResidueRMSD::ProteinModelFromString),
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

      parameters.AddInitializer
      (
        "start_ensemble_filename",
        "the name of the file containing the list of pdbs for the starting ensemble",
        util::OwnPtr< io::SerializationInterface>
        (
          new io::SerializationViaStaticFunctions< assemble::ProteinEnsemble>
          (
            command::ParameterCheckFileExistence(),
            ( &AnalyzeAccessibilityChange::EnsembleAsFilename),
            ( &AnalyzeAccessibilityChange::EnsembleFromFilename),
            &m_StartEnsemble
          )
        ),
        "start_ensemble_pdbs.ls"
      );

      parameters.AddInitializer
      (
        "normalize_by_background",
        "The background is the per residue rmsd within the ending ensemble i.e. the one provided"
        " by the main ensemble flag. If this flag is set (referring to this flag, not the main ensemble flag)"
        " then the per residue rmsd of the ending ensemble will be calculated after superimposing each model"
        " onto the provided template pdb. This per residue rmsd will then be subtracted from the per residue RMSD"
        " calculated between the two provided ensembles. This is helpful for seeing differences between the "
        " ensembles that are larger than the difference within the ending ensemble. 1=true;0=false",
        io::Serialization::GetAgent( &m_NormalizeByInternalRMSD),
        "0"
      );

      return parameters;
    }

  } // namespace restraint
} // namespace bcl

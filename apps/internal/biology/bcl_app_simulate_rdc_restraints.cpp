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
#include "app/bcl_app_apps.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_command.h"
#include "command/bcl_command_flag_dynamic.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter_check_extension.h"
#include "coord/bcl_coord_move_rotate_random.h"
#include "linal/bcl_linal_matrix_operations.h"
#include "nmr/bcl_nmr_star_rdc_handler.h"
#include "pdb/bcl_pdb_factory.h"
#include "restraint/bcl_restraint_locator_coordinates_hydrogen.h"

namespace bcl
{
  namespace app
  {

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class RestraintSimulateRdcs
    //! @brief creates an RDCDistanceAssigned restraint file.
    //!
    //! @author akinlr, weinerbe
    //!
    //! @date 08/04/2010
    //!
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API RestraintSimulateRdcs :
      public Interface
    {

    private:

    //////////
    // data //
    //////////

      //! flag for pdb file from which the restraints will be simulated
      util::ShPtr< command::FlagInterface> m_PDBFilenameFlag;

      //! flag for the specifiying an exact number of desired restraints
      util::ShPtr< command::FlagInterface> m_NumberDesiredRestraintsFlag;

      //! flag for specifying the number of restraints as a fraction of residues in the protein
      util::ShPtr< command::FlagInterface> m_FractionDesiredRestraintsFlag;

      //! flag for defining tensor values
      util::ShPtr< command::FlagStatic> m_SetTensorValuesFlag;
      util::ShPtr< command::ParameterInterface> m_TensorValueXParam;
      util::ShPtr< command::ParameterInterface> m_TensorValueYParam;

      //! flag for specifying whether or not a tensor should rotate
      util::ShPtr< command::FlagInterface> m_SetTensorRotationFlag;

      //! flag for specifying which atoms should be used in the restraints
      util::ShPtr< command::FlagInterface> m_RestraintAtomTypeFlag;

      //! flag for specifying the output file path and name
      util::ShPtr< command::FlagInterface> m_OutputFileFlag;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      RestraintSimulateRdcs();

    public:

      //! @brief Clone function
      //! @return pointer to new Quality
      RestraintSimulateRdcs *Clone() const
      {
        return new RestraintSimulateRdcs( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      // This function is used to retrieve prior names for applications that are still valid (for backwards compatibility)
      // but which should not be displayed, e.g. for help
      storage::Vector< std::string> GetDeprecatedAppNames() const
      {
        return storage::Vector< std::string>( size_t( 1), "SimulateRdcRestraints");
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief initializes the command object for that executable
      util::ShPtr< command::Command> InitializeCommand() const;

      //! @brief the Main function
      //! @return error code - 0 for success
      int Main() const;

      //! @brief takes a file and returns a vector of atom types that will be used to generate RDC data
      //! @param ATOM_TYPE_FILE a file with atom types and sequence separations
      //! @return a vector which holds all of the information in the file
      storage::Vector
      <
        storage::Pair< storage::Pair< biol::AtomType, size_t>, storage::Pair< biol::AtomType, size_t> >
      >
      GetRestraintAtomType( const std::string &ATOM_TYPE_FILE) const;

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM)
      {
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        return OSTREAM;
      }

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief creates a shptr to the RDC restraints
      //! @param PROTEIN_MODEL where list of residues and restraint distances come from
      //! @param TENSOR pass the tensor so it can be defined by the user
      //! @param ATOM_TYPE this shows what atoms the RDC should be between.  This will only work with CA and N RDCs
      //! @param NUMBER_DESIRED_RESTRAINTS how many restraints should be generated
      //! @return a ShPtr to the RDC Restraints
      static util::ShPtr< restraint::RDC> CreateRDCRestraints
      (
        const assemble::ProteinModel &PROTEIN_MODEL,
        const linal::Matrix3x3< double> &TENSOR,
        const storage::Vector
          <
            storage::Pair< storage::Pair< biol::AtomType, size_t>, storage::Pair< biol::AtomType, size_t> >
          > &ATOM_TYPES,
        const size_t NUMBER_DESIRED_RESTRAINTS
      );

    private:

      static const ApplicationType RestraintSimulateRdcs_Instance;

    }; // RestraintSimulateRdcs

    //! @brief initializes the command object for that executable
    util::ShPtr< command::Command> RestraintSimulateRdcs::InitializeCommand() const
    {
      util::ShPtr< command::Command> sp_cmd( new command::Command());

      // add flag to give PDB File
      sp_cmd->AddFlag( m_PDBFilenameFlag);

      // add flag to give the number of restraints to be given
      sp_cmd->AddFlag( m_NumberDesiredRestraintsFlag);

      // add flag that will allow user to define tensor value
      sp_cmd->AddFlag( m_SetTensorValuesFlag);

      // add flag that puts a rotation factor on the tensor
      sp_cmd->AddFlag( m_SetTensorRotationFlag);

      // add flag
      sp_cmd->AddFlag( m_RestraintAtomTypeFlag);

      // add flag
      sp_cmd->AddFlag( m_OutputFileFlag);

      // pdb factory flags
      sp_cmd->AddFlag( pdb::Factory::GetFlagMinSSESize());
      sp_cmd->AddFlag( pdb::Factory::GetFlagSSEsFromBackBone());
      pdb::Factory::GetFlagAAClass()->GetParameterList()( 0)->SetDefaultParameter( biol::GetAAClasses().e_AAComplete);
      sp_cmd->AddFlag( pdb::Factory::GetFlagAAClass());

      // add default bcl parameters
      command::GetAppDefaultFlags().AddDefaultCommandlineFlags( *sp_cmd);

      // return assembled Command object
      return sp_cmd;
    }

    //! @brief the Main function
    //! @return error code - 0 for success
    int RestraintSimulateRdcs::Main() const
    {
      // define the tensor values from the values given in the Parameter
      const double x_value( m_TensorValueXParam->GetNumericalValue< double>());
      const double y_value( m_TensorValueYParam->GetNumericalValue< double>());

      // store values in tensor vector
      linal::Matrix3x3< double> tensor;
      tensor( 0, 0) = x_value;
      tensor( 1, 1) = y_value;
      tensor( 2, 2) = -x_value - y_value;

      // apply random rotation if requested
      if( m_SetTensorRotationFlag->GetFlag())
      {
        const linal::Matrix3x3< double> rotation
        (
          coord::MoveRotateRandom::GenerateRandomRotation( linal::Vector3D( 2.0 * math::g_Pi)).GetMatrix()
        );

        tensor = rotation.Transposed() * tensor * rotation;
      }

      // open stream to write the restraints
      io::OFStream write;
      io::File::MustOpenOFStream( write, m_OutputFileFlag->GetFirstParameter()->GetValue());

      nmr::StarRDCHandler handler_rdc;
      // write RDCs to file
      handler_rdc.WriteRestraints
      (
        write,
        *CreateRDCRestraints
        (
          pdb::Factory().ProteinModelFromPDBFilename( m_PDBFilenameFlag->GetFirstParameter()->GetValue()),
          tensor,
          GetRestraintAtomType( m_RestraintAtomTypeFlag->GetFirstParameter()->GetValue()),
          m_NumberDesiredRestraintsFlag->GetFirstParameter()->GetNumericalValue< size_t>()
        )
      );

      // close the write stream
      io::File::CloseClearFStream( write);

      //successful end
      return 0;
    }

    //! @brief takes a file and returns a vector of atom types that will be used to generate RDC data
    //! @param ATOM_TYPE_FILE a file with atom types and sequence separations
    //! @return a vector which holds all of the information in the file
    storage::Vector
    <
      storage::Pair< storage::Pair< biol::AtomType, size_t>, storage::Pair< biol::AtomType, size_t> >
    >
    RestraintSimulateRdcs::GetRestraintAtomType( const std::string &ATOM_TYPE_FILE) const
    {
      // initialize write and read stream objects
      io::IFStream read;

      // open restraint atom type file
      io::File::MustOpenIFStream( read, ATOM_TYPE_FILE);

      // create a place to store the data
      storage::Vector< storage::Pair< storage::Pair< biol::AtomType, size_t>, storage::Pair< biol::AtomType, size_t> > >
        restraint_atom_type;

      // read through the file
      while( !read.eof())
      {
        // store the line
        std::string buffer;
        std::getline( read, buffer);

        // split the string to get the different things necessary for storing atom types
        const storage::Vector< std::string> split_line( util::SplitString( buffer, " "));

        if( split_line.GetSize() < 4)
        {
          break;
        }

        // create the pairs and their sequence distances
        const storage::Pair< biol::AtomType, size_t> pair_a
        (
          biol::AtomType( split_line( 0)), util::ConvertStringToNumericalValue< size_t>( split_line( 1))
        );
        const storage::Pair< biol::AtomType, size_t> pair_b
        (
          biol::AtomType( split_line( 2)), util::ConvertStringToNumericalValue< size_t>( split_line( 3))
        );

        // store the pairs in a final pair and push back into the vector
        const storage::Pair< storage::Pair< biol::AtomType, size_t>, storage::Pair< biol::AtomType, size_t> >
          pairs( pair_a, pair_b);
        restraint_atom_type.PushBack( pairs);
      }

      // close and clear the ISTREAM
      io::File::CloseClearFStream( read);
      return restraint_atom_type;
    }

    //! @brief creates a shptr to the RDC restraints
    //! @param PROTEIN_MODEL where list of residues and restraint distances come from
    //! @param TENSOR pass the tensor so it can be defined by the user
    //! @param ATOM_TYPE this shows what atoms the RDC should be between.  This will only work with CA and N RDCs
    //! @param NUMBER_DESIRED_RESTRAINTS how many restraints should be generated
    //! @return a ShPtr to the RDC Restraints
    util::ShPtr< restraint::RDC> RestraintSimulateRdcs::CreateRDCRestraints
    (
      const assemble::ProteinModel &PROTEIN_MODEL,
      const linal::Matrix3x3< double> &TENSOR,
      const storage::Vector
        <
          storage::Pair< storage::Pair< biol::AtomType, size_t>, storage::Pair< biol::AtomType, size_t> >
        > &ATOM_TYPES,
      const size_t NUMBER_DESIRED_RESTRAINTS
    )
    {
      // assert atom types not empty
      BCL_Assert( ATOM_TYPES.GetSize() != 0, "No atom types were given for creating an RDC")
      // get all amino acids
      const util::SiPtrVector< const biol::AABase> aa_vector( PROTEIN_MODEL.GetAminoAcids());

      // break if there are no amino acids
      BCL_Assert( !aa_vector.IsEmpty(), "No amino acids found");

      // initialize rdc restraint container
      restraint::RDC rdcs;

      // loop through the vector to determine which atom types to use
      for
      (
        storage::Vector< storage::Pair< storage::Pair< biol::AtomType, size_t>, storage::Pair< biol::AtomType, size_t> > >::const_iterator itr( ATOM_TYPES.Begin()),
          end_itr( ATOM_TYPES.End());
        itr != end_itr; ++itr
      )
      {
        for
        (
          util::SiPtrVector< const biol::AABase>::const_iterator aa_itr( aa_vector.Begin()),
            aa_itr_end( aa_vector.End());
            aa_itr != aa_itr_end; ++aa_itr
        )
        {
          // create locators
          restraint::LocatorCoordinatesHydrogen first_locator
          (
            ( *aa_itr)->GetChainID(),
            ( *aa_itr)->GetSeqID() + itr->First().Second(),
            itr->First().First()
          );
          restraint::LocatorCoordinatesHydrogen second_locator
          (
            ( *aa_itr)->GetChainID(),
            ( *aa_itr)->GetSeqID() + itr->Second().Second(),
            itr->Second().First()
          );

          // calculate the rdc value
          double rdc_value
          (
            restraint::RDC::CalculateValue
            (
              first_locator.Locate( PROTEIN_MODEL),
              second_locator.Locate( PROTEIN_MODEL),
              TENSOR
            )
          );

          // calculate the internuclear distance
          double internuclear_dist( first_locator.GetAtomType()->GetBondLength( second_locator.GetAtomType()));

          // if the rdc value is defined
          if( util::IsDefined( rdc_value))
          {
            // pushback the restraint
            rdcs.PushBack
            (
              first_locator,
              second_locator,
              internuclear_dist,
              rdc_value
            );
          }
        }

      }

      // set the number of restraints to write out
      size_t nr_restraints( rdcs.GetData().GetSize());

      // if the number of restraints requested is smaller than the total #
      if( NUMBER_DESIRED_RESTRAINTS < nr_restraints)
      {
        // shuffle the restraints to get a random sampling
        rdcs.Shuffle();

        // set the number of restraints
        nr_restraints = NUMBER_DESIRED_RESTRAINTS;
      }

      util::ShPtr< restraint::RDC> rdc_final( new restraint::RDC());

      // initialize an iterator

      storage::Vector< storage::Triplet< restraint::DataPairwise, double, double> >::const_iterator itr( rdcs.GetData().Begin());

      for( size_t count( 0); count != nr_restraints; count++)
      {
        rdc_final->PushBack
        (
          *( itr->First().First()),
          *( itr->First().Second()),
          itr->Second(),
          itr->Third()
        );
        ++itr;
      }

      return rdc_final;
    }

    // default constructor
    RestraintSimulateRdcs::RestraintSimulateRdcs() :
      m_PDBFilenameFlag
      (
        new command::FlagStatic
        (
          "pdb",
          "\tpdb file for which restraints will be created",
          command::Parameter
          (
            "pdb_filename",
            "\tfilename for which restraints will be calculated",
              command::ParameterCheckExtension( ".pdb"),
            ""
          )
        )
      ),
      m_NumberDesiredRestraintsFlag
      (
        new command::FlagStatic
        (
          "number_of_restraints", "Flag for the specifiying an exact number of desired restraints.",
           command::Parameter
          (
            "number_of_desired_restraints",
            "The number of restraints which are desired to be gotten.",
            util::Format()( std::numeric_limits< size_t>::max())
          )
        )
      ),
      m_SetTensorValuesFlag
      (
        new command::FlagStatic
        (
          "set_tensor_values",
          "Flag for specifying what the tensor should be"
        )
      ),
      m_TensorValueXParam
      (
        new command::Parameter
        (
          "x_tensor", "\tset X value of tensor", "-2.3"
        )
      ),
      m_TensorValueYParam
      (
        new command::Parameter
        (
          "y_tensor", "\tset Y value of tensor", "-17.3"
        )
      ),
      m_SetTensorRotationFlag
      (
        new command::FlagStatic
        (
          "rotate_tensor",
          "\tFlag for determining whether or not to rotate the tensor randomly"
        )
      ),
      m_RestraintAtomTypeFlag
      (
        new command::FlagDynamic
        (
          "atom_type", "which atoms 'N' or 'CA' should be used in the restraint",
           command::Parameter
          (
            "atom choice",
            "file formatted as [Atom Type 1] [Sequence Distance Between Atom Type 1 and2 ] "
            " [Atom Type 2] [Sequence Distance Between Atom Type 1 and2 ]",
            ""
          )
        )
      ),
      m_OutputFileFlag
      (
        new command::FlagStatic
        (
          "output_file",
          "Path and name of the output file which will hold the simulated restraints.",
          command::Parameter
          (
            "string",
            "\tpath and name of the outputfile",
            "rdc.star"
          )
        )
      )
    {
      // attach parameters to flags
      m_SetTensorValuesFlag->PushBack( m_TensorValueXParam);
      m_SetTensorValuesFlag->PushBack( m_TensorValueYParam);
    }

    const ApplicationType RestraintSimulateRdcs::RestraintSimulateRdcs_Instance
    (
      GetAppGroups().AddAppToGroup( new RestraintSimulateRdcs(), GetAppGroups().e_Restraint)
    );

  } // namespace app
} // namespace bcl

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
#include "fold/bcl_fold_loop_library.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_atom.h"
#include "command/bcl_command_command_state.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_serialization.h"
#include "linal/bcl_linal_matrix_inversion_gauss_jordan.h"
#include "pdb/bcl_pdb_factory.h"
#include "util/bcl_util_enumerated.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {

  //////////
  // data //
  //////////

    //! single instance of this class
    const util::SiPtr< const util::ObjectInterface> LoopLibrary::s_Instance
    (
      util::Enumerated< LoopLibrary>::AddInstance( new LoopLibrary())
    );

    //! map of loop libraries
    storage::HashMap< std::string, util::ShPtr< LoopLibrary> > LoopLibrary::s_LoopLibraries =
      storage::HashMap< std::string, util::ShPtr< LoopLibrary> >();

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    LoopLibrary::LoopLibrary() :
      m_DistanceBinWidth(),
      m_AngleBinWidth(),
      m_LoopTemplates(),
      m_LibraryFileName()
    {
    }

    //! @brief construct from loop template library
    //! @param LIBRARY_FILE_NAME path to the loop template library
    //! @param DISTANCE_BIN_WIDTH width of the bins to discretize the euclidean distance
    //! @param ANGLE_BIN_WIDTH width of the bins to discretize the Euler angles
    LoopLibrary::LoopLibrary
    (
      const std::string &LIBRARY_FILE_NAME,
      double DISTANCE_BIN_WIDTH,
      double ANGLE_BIN_WIDTH
    ) :
      m_DistanceBinWidth( DISTANCE_BIN_WIDTH),
      m_AngleBinWidth( ANGLE_BIN_WIDTH),
      m_LoopTemplates(),
      m_LibraryFileName( LIBRARY_FILE_NAME)
    {
      std::ostringstream oss;
      ReadInitializerSuccessHook( util::ObjectDataLabel(), oss);
    }

    //! @brief returns a loop library from the given file
    //! @param LIBRARY_FILE_NAME path to the loop template library
    //! @param DISTANCE_BIN_WIDTH width of the bins to discretize the euclidean distance
    //! @param ANGLE_BIN_WIDTH width of the bins to discretize the Euler angles
    util::ShPtr< LoopLibrary> LoopLibrary::CreateLoopLibrary
    (
      const std::string &LIBRARY_FILE_NAME,
      double DISTANCE_BIN_WIDTH,
      double ANGLE_BIN_WIDTH
    )
    {
      if( command::CommandState::GetGlobalCommandState().IsInStaticInitialization())
      {
        return util::ShPtr< LoopLibrary>();
      }
      // compute hash key for the given parameters
      std::ostringstream oss;
      oss << LIBRARY_FILE_NAME << "\t" << DISTANCE_BIN_WIDTH << "\t" << ANGLE_BIN_WIDTH;
      const std::string key( oss.str());

      // if a library with this key does not exist, create one
      if( s_LoopLibraries.Find( key) == s_LoopLibraries.End())
      {
        util::ShPtr< LoopLibrary> sp_library( new LoopLibrary( LIBRARY_FILE_NAME, DISTANCE_BIN_WIDTH, ANGLE_BIN_WIDTH));
        s_LoopLibraries[ key] = sp_library;
        return sp_library;
      }
      return s_LoopLibraries[ key];
    }

    //! @brief clone function
    //! @return pointer to a new LoopLibrary
    LoopLibrary *LoopLibrary::Clone() const
    {
      return new LoopLibrary( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns the name of this class
    //! @return the name of this class
    const std::string &LoopLibrary::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns the default distance bin width
    //! @return default distance bin width
    double LoopLibrary::GetDefaultDistanceBinWidth()
    {
      static double s_DefaultDistanceBinWidth( 2.0);
      return s_DefaultDistanceBinWidth;
    }

    //! @brief returns the default angle bin width
    //! @return default angle bin width
    double LoopLibrary::GetDefaultAngleBinWidth()
    {
      static double s_DefaultAngleBinWidth( 3.14159265359 / 180.0 * 35.0);
      return s_DefaultAngleBinWidth;
    }

    //! @brief returns the loop templates in the library
    //! @return the loop templates in the library
    const storage::HashMap< std::string, util::ShPtrVector< LoopParameters> > &LoopLibrary::GetTemplates() const
    {
      return m_LoopTemplates;
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &LoopLibrary::GetAlias() const
    {
      static const std::string s_name( "LoopLibrary");
      return s_name;
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer LoopLibrary::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "Stores loop conformations to close gaps.");
      serializer.AddInitializer
      (
        "distance bin width",
        "bin width for the coordinates of the translation vector",
        io::Serialization::GetAgent( &m_DistanceBinWidth),
        "2.0"
      );
      serializer.AddInitializer
      (
        "rotation bin width",
        "bin width for the rotation angles",
        io::Serialization::GetAgent( &m_AngleBinWidth),
        "0.61"
      );
      serializer.AddInitializer
      (
        "loop library",
        "path to the file containing the loop template library",
        io::Serialization::GetAgent( &m_LibraryFileName)
      );

      return serializer;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief read a loop template library from the given input stream
    //! @param ISTREAM input stream from which to read the loop template library
    void LoopLibrary::ReadLibrary( std::istream &ISTREAM)
    {
      BCL_MessageStd( "Reading loop template library.");

      // read in the number of templates in the library
      size_t num_templates( 0);
      ISTREAM >> num_templates;

      // read in each template and insert it into the hash map
      for( size_t i( 0); i < num_templates; ++i)
      {
        LoopParameters loop;
        ISTREAM >> loop;
        AddTemplate( loop);
      }
    }

    //! @brief extends the template library by recombining template to form larger templates
    //! @param map in the format <sequence length>:<number of new templates> to specify how many additional
    //! templates should be generated
    void LoopLibrary::RecombineTemplates( const storage::HashMap< size_t, size_t> &COUNTS)
    {
      // iterate over the given map and create the templates for each given length
      for( auto map_it( COUNTS.Begin()), map_it_end( COUNTS.End()); map_it != map_it_end; ++map_it)
      {
        // create the templates from the already existing template library
        const size_t &length( map_it->first);
        const size_t &count( map_it->second);
        util::ShPtrVector< LoopParameters> templates( GenerateTemplates( length, count));

        // add the newly created templates to the template library
        for( auto temp_it( templates.Begin()), temp_it_end( templates.End()); temp_it != temp_it_end; ++temp_it)
        {
          AddTemplate( **temp_it);
        }
      }
    }

    //! @brief returns suitable templates for the given loop parameters
    //! @param LOOP_PARAMS loop parameters for which to find templates
    //! @return loop templates for the given parameters
    util::ShPtrVector< LoopParameters> LoopLibrary::FindTemplates( const LoopParameters &LOOP_PARAMS) const
    {
      // find suitable templates or return an empty vector
      storage::HashMap< std::string, util::ShPtrVector< LoopParameters> >::const_iterator result_it
      (
        m_LoopTemplates.Find( ComputeKey( LOOP_PARAMS))
      );
      return result_it != m_LoopTemplates.End() ? result_it->second : util::ShPtrVector< LoopParameters>();
    }

    //! @brief returns suitable templates for the given loop length
    //! @param LENGTH sequence length of the loop
    //! @return loop templates for the given sequence length
    util::ShPtrVector< LoopParameters> LoopLibrary::FindTemplates( size_t LENGTH) const
    {
      // find suitable templates or return an empty vector
      storage::HashMap< size_t, util::ShPtrVector< LoopParameters> >::const_iterator result_it
      (
        m_LoopTemplatesLength.Find( LENGTH)
      );
      return result_it != m_LoopTemplatesLength.End() ? result_it->second : util::ShPtrVector< LoopParameters>();
    }

    //! @brief combines the two given templates to one larger template
    //! @param TEMPLATE_N the n-terminal template to be combined
    //! @param TEMPLATE_C the c-terminal template to be combined
    //! @return the template resulting from combination of the two given templates
    util::ShPtr< LoopParameters> LoopLibrary::Combine
    (
      const LoopParameters &TEMPLATE_N, const LoopParameters &TEMPLATE_C
    )
    {
      // compute the resulting sequence distance spanned by the combined template
      const size_t seq_dist( TEMPLATE_N.GetSequenceDistance() + TEMPLATE_C.GetSequenceDistance() + 1);

      // combine the dihedral angles of the smaller templates
      const storage::Vector< double> &angles_n( TEMPLATE_N.GetAngles());
      const storage::Vector< double> &angles_c( TEMPLATE_C.GetAngles());
      storage::Vector< double> angles( angles_n);
      angles.Append( angles_c);

      // construct the sequence
      biol::AASequence new_template( biol::GetAAClasses().e_AABackBone, seq_dist + 2, 'A');
      assemble::SSE new_loop( new_template, biol::GetSSTypes().STRAND);
      new_loop.SetToIdealConformationInPlace();
      biol::AASequenceFlexibility::SetPhiPsi
      (
        new_loop,
        1,
        storage::VectorND< 2, double>( 0.0, angles( 0)),
        biol::AASequenceFlexibility::SequenceDirection::e_CTerminal
      );
      int seq_id( 2);
      for( auto it( angles.Begin() + 1), it_end( angles.End() - 1); it != it_end; it += 2, ++seq_id)
      {
        storage::VectorND< 2, double> phi_psi( *it, *( it + 1));
        biol::AASequenceFlexibility::SetPhiPsi
        (
          new_loop,
          seq_id,
          phi_psi,
          biol::AASequenceFlexibility::SequenceDirection::e_CTerminal
        );
      }
      biol::AASequenceFlexibility::SetPhiPsi
      (
        new_loop,
        seq_dist + 2,
        storage::VectorND< 2, double>( angles( angles.GetSize() - 1), 0.0),
        biol::AASequenceFlexibility::SequenceDirection::e_CTerminal
      );

      const biol::AABase &res_next( *new_loop.GetMembers()( 1));
      const biol::Atom c_next( res_next.GetAtom( biol::GetAtomTypes().N));
      for( size_t i( 1); i <= seq_dist; ++i)
        {
          const biol::AABase &res_prev( *new_loop.GetMembers()( i-1));
          const biol::AABase &res_next( *new_loop.GetMembers()( i+1));
          const biol::Atom c_prev( res_prev.GetAtom( biol::GetAtomTypes().C));
          const biol::Atom n_prev( res_next.GetAtom( biol::GetAtomTypes().N));
        }
      const biol::AABase &res_before( *new_loop.GetMembers()( 5));
      const biol::Atom n_before( res_before.GetAtom( biol::GetAtomTypes().C));

      // construct the template from the fitted sequence
      const util::ShPtr< LoopParameters> sp_new_loop
      (
        LoopParameters::Create( *new_loop.GetFirstAA(), *new_loop.GetLastAA(), angles)
      );

      return sp_new_loop;
    }

    //! @brief creates a given number of loop templates of the given length
    //! @param LENGTH sequence length of the generated loop templates
    //! @param COUNT how many loop templates to generate
    //! @return generated loop templates
    util::ShPtrVector< LoopParameters> LoopLibrary::GenerateTemplates( size_t LENGTH, size_t COUNT) const
    {
      BCL_MessageStd( "Generating additional loop templates.");

      // the algorithm only works for template lengths >= 3
      BCL_Assert( LENGTH >= 3, "template length must be 3 or greater");

      // the templates get combined at the anchor point leading to a longer effective length
      const size_t length( LENGTH - 1);

      // generate the given number of templates of the given sequence length
      util::ShPtrVector< LoopParameters> new_templates;
      for( size_t i( 0); i < COUNT; ++i)
      {
        // try different template lengths until a fitting combination is found
        bool fit( false);
        size_t iterations( 0);
        while( !fit && iterations < 50)
        {
          // set the lengths of the templates that are to be combined
          const size_t length_1( random::GetGlobalRandom().SizeT( math::Range< size_t>( 1, length - 1)));
          const size_t length_2( length - length_1);

          // get templates of the given lengths
          storage::HashMap< size_t, util::ShPtrVector< LoopParameters> >::const_iterator result_it_1
          (
            m_LoopTemplatesLength.Find( length_1)
          );
          storage::HashMap< size_t, util::ShPtrVector< LoopParameters> >::const_iterator result_it_2
          (
            m_LoopTemplatesLength.Find( length_2)
          );
          if( result_it_1 == m_LoopTemplatesLength.End() || result_it_2 == m_LoopTemplatesLength.End())
          {
            ++iterations;
            continue;
          }
          fit = true;
          const util::ShPtrVector< LoopParameters> templates_1( result_it_1->second);
          const util::ShPtrVector< LoopParameters> templates_2( result_it_2->second);

          // randomly select a template from the found ones
          const LoopParameters template_1
          (
            **random::GetGlobalRandom().Iterator( templates_1.Begin(), templates_1.End(), templates_1.GetSize())
          );
          const LoopParameters template_2
          (
            **random::GetGlobalRandom().Iterator( templates_2.Begin(), templates_2.End(), templates_2.GetSize())
          );

          // combine the two templates and add the new templates to the result
          const util::ShPtr< LoopParameters> sp_new_template( Combine( template_1, template_2));
          new_templates.PushBack( sp_new_template);
        }
      }

      return new_templates;
    }

    //! @brief estimates how many templates of different lengths are needed to close the loops for the given
    //! protein models
    //! @param MODELS models for which the loop regions shall be constructed
    //! @return estimation of the number of additional templates required
    util::ShPtr< storage::HashMap< size_t, size_t> > LoopLibrary::EstimateTemplateRequirement
    (
      const util::ShPtrVector< assemble::ProteinModel> &MODELS
    )
    {
      BCL_MessageStd( "Estimating how many additional loop templates are needed.");

      // count the loops and their respective lengths in the given set of protein models
      util::ShPtr< storage::HashMap< size_t, size_t> > sp_estimate( new storage::HashMap< size_t, size_t>());
      for( auto model_it( MODELS.Begin()), model_it_end( MODELS.End()); model_it != model_it_end; ++model_it)
      {
        // get the loops in the current model
        const util::SiPtrVector< const assemble::SSE> loops( ( **model_it).GetSSEs( biol::GetSSTypes().COIL));

        // determine loop lengths and frequency and record them in the hash map
        for( auto loop_it( loops.Begin()), loop_it_end( loops.End()); loop_it != loop_it_end; ++loop_it)
        {
          const size_t length( ( **loop_it).GetSize());
          if( sp_estimate->Count( length) == 0)
          {
            sp_estimate->Insert( storage::Pair< size_t, size_t>( length, 1));
          }
          else
          {
            ++( *sp_estimate)[ length];
          }
        }
      }

      // for frequently occurring loops and longer loop lengths, estimate that more templates are needed
      for( auto est_it( sp_estimate->Begin()), est_it_end( sp_estimate->End()); est_it != est_it_end; ++est_it)
      {
        est_it->second = est_it->second * 10000;
      }

      return sp_estimate;
    }

    //! @brief Set the members of this object from the given LABEL
    //! @param LABEL the label containing members that should be read of this class
    //! @param ERROR_STREAM stream with which to write errors
    //! @return return code indicating success or failure
    bool LoopLibrary::ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM)
    {
      if( !command::CommandState::GetGlobalCommandState().IsInStaticInitialization())
      {
        // read the loop template library file
        io::IFStream lib_file;
        io::File::MustOpenIFStream( lib_file, m_LibraryFileName);
        ReadLibrary( lib_file);
        io::File::CloseClearFStream( lib_file);
      }

      return true;
    }

  //////////////////////
  // input and output //
  //////////////////////

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief computes the hash key for the given loop parameters
    //! @param LOOP_PARAMS loop parameters for which to compute the hash key
    //! @return hash key for the given loop parameters
    std::string LoopLibrary::ComputeKey( const LoopParameters &LOOP_PARAMS) const
    {
      // get the geometrical parameters for the computation of the hash key
      const size_t seq_dist( LOOP_PARAMS.GetSequenceDistance());
      const linal::Vector3D &translation( LOOP_PARAMS.GetTranslation());
      const linal::Vector3D &euler_angles( LOOP_PARAMS.GetRotation());

      // convert the angle bin width into radians
      const double angle_bin_width( m_AngleBinWidth / 180.0 * math::g_Pi);

      // bin the geometrical parameters
      const size_t t_x( translation( 0) / m_DistanceBinWidth);
      const size_t t_y( translation( 1) / m_DistanceBinWidth);
      const size_t t_z( translation( 2) / m_DistanceBinWidth);
      const size_t r_x( euler_angles( 0) / angle_bin_width);
      const size_t r_y( euler_angles( 1) / angle_bin_width);
      const size_t r_z( euler_angles( 2) / angle_bin_width);

      // compute the hash key from the geometrical parameters
      const std::string key
      (
        ToString( seq_dist) + ToString( t_x) + ToString( t_y) + ToString( t_z) +
        ToString( r_x) + ToString( r_y) + ToString( r_z)
      );

      return key;
    }

    //! @brief adds the given loop template to the template library
    //! @param LOOP_PARAMS loop which shall be added to the library
    void LoopLibrary::AddTemplate( const LoopParameters &LOOP_PARAMS)
    {
      // compute the hash key of the loop and add the template to the library
      const std::string key( ComputeKey( LOOP_PARAMS));
      const util::ShPtr< LoopParameters> sp_loop( util::CloneToShPtr( LOOP_PARAMS));
      m_LoopTemplates[ key].PushBack( sp_loop);

      // add the loop templates by sequence length
      const size_t length( sp_loop->GetSequenceDistance());
      m_LoopTemplatesLength[ length].PushBack( sp_loop);
    }

  } // namespace fold
} // namespace bcl

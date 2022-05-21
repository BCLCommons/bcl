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
#include "mc/bcl_mc_mutate_loop_add.h"

// includes from bcl - sorted alphabetically
#include "fold/bcl_fold_protein_geometry.h"
#include "io/bcl_io_serialization.h"
#include "score/bcl_score.h"
#include "util/bcl_util_enumerated.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace mc
  {

  //////////
  // data //
  //////////

    //! single instance of this class
    const util::SiPtr< const util::ObjectInterface> MutateLoopAdd::s_Instance
    (
      util::Enumerated< math::MutateInterface< assemble::ProteinModel> >::AddInstance( new MutateLoopAdd)
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief construct from members
    //! @param LOOP_LIBRARY_FILENAME path of the loop template library
    MutateLoopAdd::MutateLoopAdd( const std::string &LOOP_LIBRARY_FILENAME) :
      m_LoopLocator(),
      m_LoopLibraryFilename( LOOP_LIBRARY_FILENAME)
    {
      std::ostringstream oss;
      ReadInitializerSuccessHook( util::ObjectDataLabel(), oss);
    }

    //! @brief copy constructor
    //! @return pointer to a new MutateLoopAdd
    MutateLoopAdd *MutateLoopAdd::Clone() const
    {
      return new MutateLoopAdd( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns the name of this class
    //! @return the name of this class
    const std::string &MutateLoopAdd::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns the default path of the loop template library
    //! @return the default path of the loop template library
    const std::string &MutateLoopAdd::GetDefaultLibraryPath()
    {
      static std::string s_default_path( "histogram/loop_library");
      return s_default_path;
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &MutateLoopAdd::GetAlias() const
    {
      static const std::string s_name( "MutateLoopAdd");
      return s_name;
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer MutateLoopAdd::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "Constructs a missing loop region in a protein model using a template library.");
      serializer.AddInitializer
      (
        "loop library",
        "path to the loop template library",
        io::Serialization::GetAgentInputFilename( &m_LoopLibraryFilename)
      );

      return serializer;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief Set the members of this object from the given LABEL
    //! @param LABEL the label containing members that should be read of this class
    //! @param ERROR_STREAM stream with which to write errors
    //! @return return code indicating success or failure
    bool MutateLoopAdd::ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM)
    {
      m_LoopLibrary = fold::LoopLibrary::CreateLoopLibrary( m_LoopLibraryFilename);
      return true;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief applies a mutation to the given protein model
    //! @param MODEL protein model to which to apply the mutation
    //! @return result mutating the given protein model
    math::MutateResult< assemble::ProteinModel> MutateLoopAdd::operator()
    (
      const assemble::ProteinModel &MODEL
    ) const
    {
      // model resulting from applying the mutate
      util::ShPtr< assemble::ProteinModel> result_model( MODEL.Clone());

      // find missing loop regions in the protein model
      const util::ShPtrVector< fold::LoopParameters> loops( m_LoopLocator.Locate( MODEL));

      // if there are no missing loops in the protein model return the original model
      if( loops.IsEmpty())
      {
        BCL_MessageVrb( "No missing loops in the protein model. Returning original model.");
        return math::MutateResult< assemble::ProteinModel>( result_model, *this);
      }

      // randomly select one of the missing loops for construction
      const size_t loop_index( random::GetGlobalRandom().SizeT( math::Range< size_t>( 0, loops.GetSize() - 1)));
      const fold::LoopParameters &loop( *loops( loop_index));

      // find a suitable loop template in the library
      const util::ShPtrVector< fold::LoopParameters> templates( m_LoopLibrary->FindTemplates( loop));

      // loop region can't be constructed if there is no suitable template
      if( templates.IsEmpty())
      {
        BCL_MessageVrb( "No suitable loop template found. Returning original model.");
        return math::MutateResult< assemble::ProteinModel>( result_model, *this);
      }

      // select a random template from the found templates
      const util::ShPtrVector< fold::LoopParameters>::const_iterator selected_template_it
      (
        random::GetGlobalRandom().Iterator( templates.Begin(), templates.End(), templates.GetSize())
      );
      const fold::LoopParameters &selected_template( **selected_template_it);

      // fit the loop to the selected template and add it to the model
      result_model = fold::ProteinGeometry::FitToTemplate( MODEL, loop, selected_template);

      // return the result
      return math::MutateResult< assemble::ProteinModel>( result_model, *this);
    }

  //////////////////////
  // input and output //
  //////////////////////

  } // namespace mc
} // namespace bcl

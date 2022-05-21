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
#include "assemble/bcl_assemble_protein_ensemble.h"
#include "fold/bcl_fold_mutate_protein_model_replicate_conformation.h"
#include "math/bcl_math_mutate_result.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> MutateProteinModelReplicateConformation::s_Instance
    (
      GetObjectInstances().AddInstance( new MutateProteinModelReplicateConformation())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MutateProteinModelReplicateConformation::MutateProteinModelReplicateConformation() :
      m_ConformationCollector(),
      m_NumberReplications(),
      m_Scheme( GetStaticClassName< MutateProteinModelReplicateConformation>())
    {
    }

    //! @brief constructor taking parameters
    //! @param LOCATOR the method for locating the conformation to replicate
    //! @param NUMBER_REPLICATIONS the number of times each collected conformation should be copied
    //! @param SCHEME Scheme to be used
    MutateProteinModelReplicateConformation::MutateProteinModelReplicateConformation
    (
      const util::ShPtr
      <
        find::CollectorInterface< util::SiPtrList< const assemble::ProteinModel>, assemble::ProteinModel>
      > &CONFORMATION,
      int NUMBER_REPLICATIONS,
      const std::string &SCHEME
    ) :
      m_ConformationCollector( CONFORMATION),
      m_NumberReplications( NUMBER_REPLICATIONS),
      m_Scheme( SCHEME)
    {
    }

    //! @brief Clone function
    //! @return pointer to new MutateProteinModelReplicateConformation
    MutateProteinModelReplicateConformation *MutateProteinModelReplicateConformation::Clone() const
    {
      return new MutateProteinModelReplicateConformation( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &MutateProteinModelReplicateConformation::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////
  
  ///////////////
  // operators //
  ///////////////

    //! @brief virtual operator taking an ARGUMENT and returning a mutate object of t_ArgumentType
    //! @param PROTEIN_MODEL protein model interest
    //! @return MutateResult with ProteinModel after the mutate
    math::MutateResult< assemble::ProteinModel>
    MutateProteinModelReplicateConformation::operator()( const assemble::ProteinModel &PROTEIN_MODEL) const
    {
      // static empty model
      const math::MutateResult< assemble::ProteinModel> empty_result( util::ShPtr< assemble::ProteinModel>(), *this);

      // try to locate the desired models
      util::SiPtrList< const assemble::ProteinModel> conformations( m_ConformationCollector->Collect( PROTEIN_MODEL));

      // if the desired conformation could not be found
      if( conformations.IsEmpty())
      {
        return empty_result;
      }

      // create new model from the first element located conformation
      util::ShPtr< assemble::ProteinEnsemble> new_model
      (
        new assemble::ProteinEnsemble( *conformations.FirstElement())
      );

      // copy the first model one less than the desired number of times, since it is already set as the model
      for( int copy( 0); copy < m_NumberReplications; ++copy)
      {
        // make new copy of the conformation
        util::ShPtr< assemble::ProteinModel> model_copy( conformations.FirstElement()->HardCopy());

        new_model->InsertElement( model_copy);
      }

      // iterate through the rest of the conformations
      for
      (
        util::SiPtrList< const assemble::ProteinModel>::const_iterator
          model_itr( ++conformations.Begin()), //< start at the second element
          model_itr_end( conformations.End());
        model_itr != model_itr_end;
        ++model_itr
      )
      {
        // copy the rest of the models the desired number of times
        for( int copy( 0); copy < m_NumberReplications + 1; ++copy)
        {
          new_model->InsertElement( util::ShPtr< assemble::ProteinModel>( ( *model_itr)->HardCopy()));
        }
      }

      return math::MutateResult< assemble::ProteinModel>( new_model, *this);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &MutateProteinModelReplicateConformation::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_ConformationCollector, ISTREAM);
      io::Serialize::Read( m_NumberReplications, ISTREAM);
      io::Serialize::Read( m_Scheme, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &MutateProteinModelReplicateConformation::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_ConformationCollector, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_NumberReplications, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Scheme, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

  } // namespace fold
} // namespace bcl

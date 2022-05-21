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
#include "fold/bcl_fold_mutate_protein_model_filter_conformations.h"
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
    const util::SiPtr< const util::ObjectInterface> MutateProteinModelFilterConformations::s_Instance
    (
      GetObjectInstances().AddInstance( new MutateProteinModelFilterConformations())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MutateProteinModelFilterConformations::MutateProteinModelFilterConformations() :
      m_ConformationCollector(),
      m_Scheme( GetStaticClassName< MutateProteinModelFilterConformations>())
    {
    }

    //! @brief constructor taking parameters
    //! @param LOCATOR the method for locating an alternative conformation
    //! @param SCHEME Scheme to be used
    MutateProteinModelFilterConformations::MutateProteinModelFilterConformations
    (
      const util::ShPtr
      <
        find::CollectorInterface< util::SiPtrList< const assemble::ProteinModel>, assemble::ProteinModel>
      > &COLLECTOR,
      const std::string &SCHEME
    ) :
      m_ConformationCollector( COLLECTOR),
      m_Scheme( SCHEME)
    {
    }

    //! @brief Clone function
    //! @return pointer to new MutateProteinModelFilterConformations
    MutateProteinModelFilterConformations *MutateProteinModelFilterConformations::Clone() const
    {
      return new MutateProteinModelFilterConformations( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &MutateProteinModelFilterConformations::GetClassIdentifier() const
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
    MutateProteinModelFilterConformations::operator()( const assemble::ProteinModel &PROTEIN_MODEL) const
    {
      // static empty model
      const math::MutateResult< assemble::ProteinModel> empty_result( util::ShPtr< assemble::ProteinModel>(), *this);

      // collect the desired models
      const util::SiPtrList< const assemble::ProteinModel> filtered_models
      (
        m_ConformationCollector->Collect( PROTEIN_MODEL)
      );

      // if no models were collected
      if( filtered_models.IsEmpty())
      {
        return empty_result;
      }

      // initialize new ensemble using the first collected model
      util::ShPtr< assemble::ProteinEnsemble> model( new assemble::ProteinEnsemble( *filtered_models.FirstElement()));

      // insert the rest of the collected models
      for
      (
        util::SiPtrList< const assemble::ProteinModel>::const_iterator
          model_itr( ++filtered_models.Begin()), //< start at the second element
          model_itr_end( filtered_models.End());
        model_itr != model_itr_end;
        ++model_itr
      )
      {
        model->InsertElement( util::ShPtr< assemble::ProteinEnsemble>( ( *model_itr)->HardCopy()));
      }

      // return the mutated model
      return math::MutateResult< assemble::ProteinModel>( model, *this);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &MutateProteinModelFilterConformations::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_ConformationCollector, ISTREAM);
      io::Serialize::Read( m_Scheme, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &MutateProteinModelFilterConformations::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_ConformationCollector, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Scheme, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

  } // namespace fold
} // namespace bcl

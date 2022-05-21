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
#include "chemistry/bcl_chemistry_score_function_generic.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_fragment_complete.h"
#include "io/bcl_io_serialization.h"
#include "io/bcl_io_serializer.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> ScoreFunctionGeneric::s_Instance
    (
      GetObjectInstances().AddInstance( new ScoreFunctionGeneric())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    ScoreFunctionGeneric::ScoreFunctionGeneric() :
        m_Descriptor()
    {
    }

    //! @brief constructor with parameters
    //! @param DESCRIPTOR the descriptor to use
    ScoreFunctionGeneric::ScoreFunctionGeneric
    (
      const descriptor::CheminfoProperty &DESCRIPTOR
    ) :
      m_Descriptor( DESCRIPTOR)
    {
    }

    //! @brief Clone function
    //! @return pointer to new ScoreFunctionGeneric
    ScoreFunctionGeneric *ScoreFunctionGeneric::Clone() const
    {
      return new ScoreFunctionGeneric( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ScoreFunctionGeneric::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get the class name when used in a dynamic context
    //! @return the class name when used in a dynamic context
    const std::string &ScoreFunctionGeneric::GetAlias() const
    {
      static const std::string s_alias( "ScoreFunctionGeneric");
      return s_alias;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief evaluate clashes for given atom pair
    //! @param MOLECULE molecule that needs to scored
    //! @return score of MOLECULE
    double ScoreFunctionGeneric::operator()
    (
      const FragmentComplete &MOLECULE
    ) const
    {
      // initialize activity
      double activity( util::GetUndefinedDouble());

      // setup score function options
      if( m_Descriptor.IsDefined())
      {
        // use passed property
        linal::Vector< double> properties( m_Descriptor->SumOverObject( MOLECULE));
        activity = properties.Sum() / properties.GetSize();
      }
      else
      {
        // flat score landscape
        BCL_MessageStd( "No score defined; returning 0.0 to approximator!");
        activity = 0.0;
      }

      //end
      return activity;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief gets a serializer for constructing this class in a dynamic context
    //! @return the serializer containing member data
    io::Serializer ScoreFunctionGeneric::GetSerializer() const
    {
      io::Serializer member_data;

      member_data.SetClassDescription( "scores molecules using the raw mean output from a descriptor");

      member_data.AddInitializer
      (
        "descriptor",
        "the descriptor to calculate; if multi-valued, this will return the mean value.",
        io::Serialization::GetAgent( &m_Descriptor)
      );
      return member_data;
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ScoreFunctionGeneric::Read( std::istream &ISTREAM)
    {
      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &ScoreFunctionGeneric::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      return OSTREAM;
    }

  } // namespace chemistry
} // namespace bcl

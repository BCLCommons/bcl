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

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "chemistry/bcl_chemistry_pick_fragment_property_weighted.h"
#include "util/bcl_util_string_functions.h"

namespace bcl
{
  namespace chemistry
  {
  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> PickFragmentPropertyWeighted::s_Instance
    (
      GetObjectInstances().AddInstance( new PickFragmentPropertyWeighted())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! default constructor
    PickFragmentPropertyWeighted::PickFragmentPropertyWeighted() :
       m_PropertyName()
    {
    }

    //! @brief constructor from property name
    //! @param PROPERTY the name of the property over which to weight
    PickFragmentPropertyWeighted::PickFragmentPropertyWeighted( const std::string &PROPERTY_NAME) :
      m_PropertyName( PROPERTY_NAME)
    {
    }

    //! @brief Clone is the virtual Clone constructor
    //! @return a pointer to new LocatorAtomFurthest which is a copy of this
    PickFragmentPropertyWeighted *PickFragmentPropertyWeighted::Clone() const
    {
      return new PickFragmentPropertyWeighted( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &PickFragmentPropertyWeighted::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! Picks a random fragment
    //! @param FRAGMENTS is the list of fragments
    //! @return returns a random fragment
    const FragmentComplete &PickFragmentPropertyWeighted::Pick( const FragmentEnsemble &FRAGMENTS) const
    {
      BCL_Assert( FRAGMENTS.GetSize() > 0, "Fragment list is empty!");

      // cumulative_weights is the cumulative sum of the property
      storage::List< float> cumulative_weights;
      float cumulative_weight( 0.0);

      for
      (
        FragmentEnsemble::const_iterator itr( FRAGMENTS.Begin()), itr_end( FRAGMENTS.End());
        itr != itr_end;
        ++itr
      )
      {
        float weight( 1.0);
        const std::string &property( itr->GetMDLProperty( m_PropertyName));
        // empty string if property is missing

        if( !property.size())
        {
          BCL_MessageStd( "Property " + m_PropertyName + " not found for fragment.");
        }
        else if( !util::TryConvertFromString( weight, util::TrimString( property), util::GetLogger()))
        {
          BCL_MessageStd( "Property " + m_PropertyName + " has value not convertible to number: '" + property + "'");
        }

        cumulative_weight += math::Absolute( weight);
        cumulative_weights.PushBack( cumulative_weight);
      }

      float threshold( random::GetGlobalRandom().Random( float( 0), cumulative_weight));

      // go through the cumulative_weights list until we exceed the threshold weight
      FragmentEnsemble::const_iterator frag_itr( FRAGMENTS.Begin()), frag_itr_end( FRAGMENTS.End());

      for
      (
        storage::List< float>::const_iterator itr( cumulative_weights.Begin()), itr_end( cumulative_weights.End());
        itr != itr_end && frag_itr != frag_itr_end && *itr < threshold;
        ++itr, ++frag_itr
      )
      {
      }

      if( frag_itr == frag_itr_end)
      {
        // something went wrong - just pick a random fragment
        return *random::GetGlobalRandom().Iterator( FRAGMENTS.Begin(), FRAGMENTS.End(), FRAGMENTS.GetSize());
      }

      // return
      return *frag_itr;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &PickFragmentPropertyWeighted::Read( std::istream &ISTREAM)
    {
      // end
      return ISTREAM;
    }

    //! @brief read from std::ostream
    //! @param OSTREAM input stream
    //! @param INDENT number of indentations
    //! @return ostream which was read from
    std::ostream &PickFragmentPropertyWeighted::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // end
      return OSTREAM;
    }

  } // namespace chemistry
} // namespace bcl

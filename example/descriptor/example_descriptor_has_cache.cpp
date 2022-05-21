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

// include example header
#include "example.h"
// include the header of the class which this example is for
#include "descriptor/bcl_descriptor_has_cache.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_descriptor_has_cache.cpp
  //! @details Tests HasCache mixin class, and indirectly CacheMap
  //!
  //! @author mendenjl
  //! @date Nov 9, 2012
  //! @remarks status complete
  //! @remarks
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDescriptorHasCache :
    public ExampleInterface
  {
  public:

    ExampleDescriptorHasCache *Clone() const
    {
      return new ExampleDescriptorHasCache( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    int Run() const
    {

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test that the mixin class performs the same functionality
      descriptor::HasCache< util::ObjectDataLabel> trivial_cache;

    /////////////////
    // data access //
    /////////////////

      util::ObjectDataLabel example_label( "example");
      linal::Vector< float> example_values( size_t( 1), 3.14159);
      BCL_ExampleCheck( trivial_cache.IsCached( example_label), false);
      BCL_ExampleCheck( trivial_cache.GetFromCache( example_label), linal::Vector< float>());
      BCL_ExampleCheck( trivial_cache.FindInCache( example_label).IsDefined(), false);

      trivial_cache.Cache( example_label, example_values);
      BCL_ExampleIndirectCheck( trivial_cache.IsCached( example_label), true, "Cache/IsCached");
      BCL_ExampleCheck( trivial_cache.GetFromCache( example_label), example_values);
      BCL_ExampleCheck( trivial_cache.FindInCache( example_label).IsDefined(), true);

      util::ObjectDataLabel example_sorted( "example(alpha=0.1,beta=0.2)");

      // check that the whole label is compared
      BCL_ExampleIndirectCheck( trivial_cache.IsCached( example_sorted), false, "Whole label comparison");

      // check uncache
      trivial_cache.Uncache( example_label);
      BCL_ExampleIndirectCheck( trivial_cache.IsCached( example_label), false, "Uncache");

      // now cache the example_sorted label
      trivial_cache.Cache( example_sorted, example_values);

      util::ObjectDataLabel example_unsorted( "example(beta = 0.2, alpha = 0.1)");
      // test that find works even though the parameters are unsorted
      BCL_ExampleIndirectCheck
      (
        trivial_cache.IsCached( example_unsorted),
        true,
        "Cache works independently of argument order"
      );

      // check that sequence ordering is respected
      util::ObjectDataLabel example_sequence_alpha_beta( "example(alpha,beta)");
      util::ObjectDataLabel example_sequence_beta_alpha( "example(beta,alpha)");

      trivial_cache.Cache( example_sequence_alpha_beta, example_values);
      BCL_ExampleIndirectCheck
      (
        trivial_cache.IsCached( example_sequence_beta_alpha),
        false,
        "Cache respects sequence order"
      );

      // check reset cache
      trivial_cache.ResetCache();
      BCL_ExampleIndirectCheck( trivial_cache.IsCached( example_unsorted), false, "ResetCache");

      // create a map from string to strings to test CacheNumeric
      storage::Map< std::string, std::string> map_str_to_value;
      map_str_to_value[ "One"] = "1";
      map_str_to_value[ "Primes(3)"] = "2 3 5";
      map_str_to_value[ example_sorted.ToString()] = "1.5";
      map_str_to_value[ "NaN"] = "nan";
      map_str_to_value[ "Letters"] = "ABC";

      trivial_cache.CacheNumeric( map_str_to_value);

      BCL_ExampleIndirectCheck( trivial_cache.IsCached( util::ObjectDataLabel( "One")), true, "CacheNumeric");
      BCL_ExampleIndirectCheck
      (
        trivial_cache.GetFromCache( util::ObjectDataLabel( "One")),
        linal::Vector< float>( size_t( 1), 1.0),
        "CacheNumeric"
      );
      BCL_ExampleIndirectCheck( trivial_cache.IsCached( util::ObjectDataLabel( "Primes(3)")), true, "CacheNumeric");
      BCL_ExampleIndirectCheck( trivial_cache.IsCached( example_sorted), true, "CacheNumeric");
      BCL_ExampleIndirectCheck( trivial_cache.IsCached( util::ObjectDataLabel( "NaN")), true, "CacheNumeric");
      BCL_ExampleIndirectCheck( trivial_cache.IsCached( util::ObjectDataLabel( "Letters")), false, "CacheNumeric");
      BCL_ExampleIndirectCheck
      (
        trivial_cache.GetFromCache( example_sorted),
        linal::Vector< float>( size_t( 1), 1.5),
        "CacheNumeric"
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleSmallMolecule

  const ExampleClass::EnumType ExampleDescriptorHasCache::s_Instance
  (
    GetExamples().AddEnum( ExampleDescriptorHasCache())
  );

} // namespace bcl

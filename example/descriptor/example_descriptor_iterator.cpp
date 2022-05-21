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
#include "descriptor/bcl_descriptor_iterator.h"

// includes from bcl - sorted alphabetically
#include "descriptor/bcl_descriptor_example_string_sequence.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //
  //! @example example_descriptor_iterator.cpp
  //!
  //! @author mendenjl
  //! @date Dec 12, 2012
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDescriptorIterator :
    public ExampleInterface
  {
  public:

    ExampleDescriptorIterator *Clone() const
    {
      return new ExampleDescriptorIterator( *this);
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

    //! @brief get the iteration result as a string (semi-colons between entries)
    //! @param TYPE type of sequence iterator to create
    //! @param STRING the string of interest
    //! @return the resulting string
    std::string WriteIterations( const descriptor::Type &TYPE, const std::string &STRING) const
    {
      descriptor::Iterator< char> itr( TYPE);
      descriptor::StringSequence string_seq( STRING);
      itr.SetObject( string_seq);
      std::ostringstream strstrm;
      size_t n_itr( 0);
      storage::Vector< size_t> positions( TYPE.GetDimension(), size_t( 0));

      // keep track of whether TYPE.GetPosition was always consistent with iterator.GetPosition
      // the iterator tracks its position directly, while TYPE.GetPosition takes the positions of the iterator and is
      // necessary when constructing descriptor::Iterators from multiple generic iterators
      bool iterator_get_position_equals_type_get_position( true);
      for( ; itr.NotAtEnd(); ++itr, ++n_itr)
      {
        size_t result_index( 0);
        for
        (
          descriptor::Iterator< char>::const_iterator
            itr_char( itr.Begin()), itr_char_end( itr.End());
          itr_char != itr_char_end;
          ++itr_char, ++result_index
        )
        {
          positions( result_index) = itr_char->GetPosition();
          strstrm << **itr_char;
        }
        strstrm << "; ";
        if( TYPE.GetPosition( positions, STRING.size()) != itr.GetPosition())
        {
          iterator_get_position_equals_type_get_position = false;
          BCL_ExampleCheck( TYPE.GetPosition( positions, STRING.size()), itr.GetPosition());
        }
      }
      BCL_ExampleIndirectCheck
      (
        iterator_get_position_equals_type_get_position,
        true,
        "TYPE.GetPosition( positions, STRING.size()) == itr.GetPosition())"
      );
      BCL_ExampleCheck( n_itr, itr.GetSize());
      BCL_ExampleCheck( itr.GetPosition(), itr.GetSize());

      return strstrm.str();
    }

    int Run() const
    {
      // create a bunch of different types of varied complexity
      descriptor::Type scalar( 0, true, descriptor::Type::e_Asymmetric);
      descriptor::Type elementwise( 1, true, descriptor::Type::e_Asymmetric);
      descriptor::Type pairs_all( 2, true, descriptor::Type::e_Asymmetric);
      descriptor::Type triplets_all( 3, true, descriptor::Type::e_Asymmetric);
      descriptor::Type triplets_no_repeat( 3, false, descriptor::Type::e_Asymmetric);
      descriptor::Type pair_all_combos( 2, true, descriptor::Type::e_Symmetric);

      // create a test string
      std::string string_to_test( "ABCD");

      // check that the correct sequence tuples are visited for each descriptor type
      BCL_ExampleCheck( WriteIterations( scalar, string_to_test), "; ");
      BCL_ExampleCheck( WriteIterations( elementwise, string_to_test), "A; B; C; D; ");
      BCL_ExampleCheck
      (
        WriteIterations( pairs_all, string_to_test),
        "AA; AB; AC; AD; BA; BB; BC; BD; CA; CB; CC; CD; DA; DB; DC; DD; "
      );

      BCL_ExampleCheck
      (
        WriteIterations( triplets_all, string_to_test),
        "AAA; AAB; AAC; AAD; ABA; ABB; ABC; ABD; ACA; ACB; ACC; ACD; ADA; ADB; ADC; ADD; "
        "BAA; BAB; BAC; BAD; BBA; BBB; BBC; BBD; BCA; BCB; BCC; BCD; BDA; BDB; BDC; BDD; "
        "CAA; CAB; CAC; CAD; CBA; CBB; CBC; CBD; CCA; CCB; CCC; CCD; CDA; CDB; CDC; CDD; "
        "DAA; DAB; DAC; DAD; DBA; DBB; DBC; DBD; DCA; DCB; DCC; DCD; DDA; DDB; DDC; DDD; "
      );

      BCL_ExampleCheck
      (
        WriteIterations( triplets_no_repeat, string_to_test),
        "ABC; ABD; ACB; ACD; ADB; ADC; "
        "BAC; BAD; BCA; BCD; BDA; BDC; "
        "CAB; CAD; CBA; CBD; CDA; CDB; "
        "DAB; DAC; DBA; DBC; DCA; DCB; "
      );

      BCL_ExampleCheck
      (
        WriteIterations( pair_all_combos, string_to_test),
        "AA; AB; AC; AD; BB; BC; BD; CC; CD; DD; "
      );

      descriptor::Type triplet_all_combos( 3, true, descriptor::Type::e_Symmetric);
      BCL_ExampleCheck
      (
        WriteIterations( triplet_all_combos, string_to_test),
        "AAA; AAB; AAC; AAD; ABB; ABC; ABD; ACC; ACD; ADD; "
        "BBB; BBC; BBD; BCC; BCD; BDD; CCC; CCD; CDD; DDD; "
      );

      descriptor::Type quartet_all_combos( 4, true, descriptor::Type::e_Symmetric);
      BCL_ExampleCheck
      (
        WriteIterations( quartet_all_combos, string_to_test),
        "AAAA; AAAB; AAAC; AAAD; AABB; AABC; AABD; AACC; AACD; AADD; "
        "ABBB; ABBC; ABBD; ABCC; ABCD; ABDD; ACCC; ACCD; ACDD; ADDD; "
        "BBBB; BBBC; BBBD; BBCC; BBCD; BBDD; BCCC; BCCD; BCDD; BDDD; "
        "CCCC; CCCD; CCDD; CDDD; DDDD; "
      );

      descriptor::Type quartet_unique_combos( 4, false, descriptor::Type::e_Symmetric);
      BCL_ExampleCheck
      (
        WriteIterations( quartet_unique_combos, string_to_test),
        "ABCD; "
      );

      descriptor::Type quartet_unique_perms( 4, false, descriptor::Type::e_Asymmetric);
      BCL_ExampleCheck
      (
        WriteIterations( quartet_unique_perms, string_to_test),
        "ABCD; ABDC; ACBD; ACDB; ADBC; ADCB; "
        "BACD; BADC; BCAD; BCDA; BDAC; BDCA; "
        "CABD; CADB; CBAD; CBDA; CDAB; CDBA; "
        "DABC; DACB; DBAC; DBCA; DCAB; DCBA; "
      );

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleDescriptorIterator

  const ExampleClass::EnumType ExampleDescriptorIterator::s_Instance
  (
    GetExamples().AddEnum( ExampleDescriptorIterator())
  );

} // namespace bcl

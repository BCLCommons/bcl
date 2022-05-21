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

#ifndef BCL_DESCRIPTOR_TYPE_H_
#define BCL_DESCRIPTOR_TYPE_H_

// include the namespace header
#include "bcl_descriptor.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "storage/bcl_storage_vector.h"
#include "util/bcl_util_object_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Type
    //! @brief information about the type of a sequence or scalar descriptor
    //! State info is currently restricted to information about descriptors of complete object (1 training example per object)
    //! or 1 training example per N-tuple of objects
    //!
    //! @see @link example_descriptor_type.cpp @endlink
    //! @author mendenjl
    //! @date Nov 26, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API Type :
      public util::ObjectInterface
    {

    public:

    //////////
    // enum //
    //////////

      //! type of symmetry (for descriptors with dimension > 1)
      //! This is an enum because additional types of symmetry are available in 3+ dimensional permutations,
      //! @see @link http://home.roadrunner.com/~hinnant/combinations.html @endlink
      //! gives 2 additional types (Circular, Reversible) in 3 dimensions
      //! In 4+ dimensions, one can also have circular reversible (this is equivalent to symmetric in dimension 3)
      //! as well as other symmetry groups
      //! For simplicity, however, only types that are currently needed are shown here
      enum Symmetry
      {
        e_Asymmetric,      //!< Order of sub-objects matter, e.g. A-B is not the same as B-A
        e_Symmetric,       //!< Order of sub-objects is irrelevant, e.g. A+B is identical to B+A
        s_NumberSymmetries //!< Number of different symmetry types
      };

    //////////
    // data //
    //////////

    private:

      size_t   m_Dimension;       //!< Number of objects from within the sequence to consider
      bool     m_ConsiderRepeats; //!< Consider descriptors with repeated elements, e.g. B,A,B
      Symmetry m_Symmetry;        //!< Symmetry, see comment above

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      Type();

      //! @brief constructor from members: dimension, whether to consider repeats, and symmetry
      Type( const size_t DIMENSION, const bool &REPEATS, const Symmetry &SYMMETRY);

      //! @brief Clone function
      //! @return pointer to new Type
      Type *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief Get the symmetry of the descriptor
      //! @return the symmetry of the descriptor
      const Symmetry &GetSymmetry() const;

      //! @brief get whether to consider repeated elements for training examples
      //! @detail for example, whether there should be a feature for elements B, B
      //! @return true if there should be features for repeated elements
      const bool &ConsiderRepeatedObjects() const;

      //! @brief get the dimension of the descriptor
      //! @return the # of sub objects that go into calculating the descriptor
      const size_t &GetDimension() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief get the number of features for the object, given the number of sub objects
      //! @param SIZE the number of sub objects for the object (e.g. atoms in the molecule, AAs in the sequences)
      //! @return the number of features for the object given the dimension and symmetry
      size_t GetNumberFeatures( const size_t &SIZE) const;

      //! @brief given a vector of element positions, return the overall position for the iterator of this type
      //! @param POSITIONS vector of size_t positions for each element
      //! @param SIZE the number of sub objects for the object (e.g. atoms in the molecule, AAs in the sequences)
      //! @return the overall position for the iterator of this type
      size_t GetPosition( const storage::Vector< size_t> &POSITIONS, const size_t &SIZE) const;

      //! @brief generalize this type to allow it to generate descriptors appropriate for another type
      //! @details for example become asymmetric if the other descriptor is asymmetric
      //! @param OTHER the type to generalize for
      //! @return reference to this
      Type &GeneralizeToHandle( const Type &OTHER);

      //! @brief equality test
      //! @param OTHER other type to test for equality
      bool operator ==( const Type &OTHER) const;

      //! @brief inequality test
      //! @param OTHER other type to test for inequality
      bool operator !=( const Type &OTHER) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class Type

  } // namespace descriptor
} // namespace bcl

#endif // BCL_DESCRIPTOR_TYPE_H_

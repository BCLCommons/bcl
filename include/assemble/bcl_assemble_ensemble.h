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

#ifndef BCL_ASSEMBLE_ENSEMBLE_H_
#define BCL_ASSEMBLE_ENSEMBLE_H_

// include the namespace header
#include "bcl_assemble.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"
#include "util/bcl_util_serializable_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Ensemble
    //! @brief Class for storing data ensembles
    //! @detail Stores a data ensemble with associated information.
    //!
    //! @see @link example_assemble_ensemble.cpp @endlink
    //! @author fischea
    //! @date Oct 29, 2016
    //!
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////
    template< typename t_DataType>
    class Ensemble :
      public util::SerializableInterface
    {

    ///////////
    // types //
    ///////////

    public:

      /////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //!
      //! @class Element
      //! @brief Stores one element of the ensemble alongside corresponding information
      //!
      //! @remarks example unnecessary
      //! @author fischea
      //! @date Nov 5, 2016
      //!
      /////////////////////////////////////////////////////////////////////////////////////////////////////////////
      class Element :
        public storage::Pair< t_DataType, double>,
        public util::SerializableInterface
      {

      private:

        typedef storage::Pair< t_DataType, double> Base;

      public:

        //! @brief default constructor
        Element() :
          Base()
        {
        }

        //! @brief construct from ensemble element and corresponding information
        //! @param ELEMENT this element of the ensemble
        //! @param POPULATION_SIZE population size of this element
        Element( t_DataType &ELEMENT, double POPULATION_SIZE = 1.0) :
          Base( ELEMENT, POPULATION_SIZE)
        {
        }

        //! @brief clone function
        //! @return pointer to a new Element
        Element *Clone() const
        {
          return new Element( *this);
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

        //! @brief return the ensemble element
        //! @return the ensemble element
        t_DataType &GetElement()
        {
          return this->First();
        }

        //! @brief return the element's population size
        //! @return the element's population size
        double GetPopulationSize() const
        {
          return this->Second();
        }

        //! @brief alias
        //! @return an alias for the class
        const std::string &GetAlias() const
        {
          static const std::string s_empty;
          return s_empty;
        }

        //! @brief return parameters for member data that are set up from the labels
        //! @return parameters for member data that are set up from the labels
        io::Serializer GetSerializer() const
        {
          io::Serializer parameters;
          parameters.SetClassDescription( "ensemble element and corresponding information");

          return parameters;
        }

      }; // class Element

    //////////////
    // typedefs //
    //////////////

      //! iterator for the members of the ensemble
      typedef typename storage::Vector< Element>::iterator iterator;

    //////////
    // data //
    //////////

    public:

      //! single instance of this class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    private:

      //! elements of the ensemble
      storage::Vector< Element> m_Elements;

    ////////////////////////////////////
    //  construction and destruction  //
    ////////////////////////////////////

    public:

      //! @brief default constructor
      Ensemble() :
        m_Elements()
      {
      }

      //! @brief construct from element iterators
      //! @param ITERATOR iterator pointing to the first element of the ensemble
      //! @param ITERATOR_END iterator pointing to the end of the ensemble
      template< typename t_InputIterator>
      Ensemble( const t_InputIterator &ITERATOR, const t_InputIterator &ITERATOR_END) :
        m_Elements( ITERATOR, ITERATOR_END)
      {
      }

      //! @brief clone function
      //! @return pointer to a new Ensemble
      Ensemble *Clone() const
      {
        return new Ensemble( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns the name of this class
      //! @return the name of this class
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief get the name of the object when used in a dynamic context
      //! @return the name of the object when used in a dynamic context
      const std::string &GetAlias() const
      {
        static const std::string s_alias( "Ensemble");
        return s_alias;
      }

      //! @brief returns parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const
      {
        io::Serializer serializer;
        serializer.SetClassDescription( "Ensemble of data instances.");

        return serializer;
      }

      //! @brief returns the size of the ensemble
      //! @return the number of elements in the ensemble
      size_t GetSize() const
      {
        return m_Elements.GetSize();
      }

      //! @brief return iterator to the first member of the ensemble
      //! @return iterator to the first member of the ensemble
      iterator Begin()
      {
        return m_Elements.Begin();
      }

      //! @brief return iterator to the end of the ensemble (i.e. behind the last element)
      //! @return iterator to the end of the ensemble
      iterator End()
      {
        return m_Elements.End();
      }

      //! @brief remove one element from the ensemble
      //! @param ITERATOR iterator to the element to be removed
      void RemoveElement( iterator ITERATOR)
      {
        m_Elements.RemoveElement( ITERATOR);
      }

      //! @brief add element to the ensemble
      //! @param ELEMENT element to add to the ensemble
      //! @param POPULATION_SIZE population size of the element
      void AddElement( t_DataType &ELEMENT, double POPULATION_SIZE = 1.0)
      {
        Element element( ELEMENT, POPULATION_SIZE);
        m_Elements.PushBack( element);
      }

      //! @brief add element to the ensemble
      //! @param ELEMENT element to add to the ensemble
      void AddElement( Element &ELEMENT)
      {
        m_Elements.PushBack( ELEMENT);
      }

    ///////////////
    // operators //
    ///////////////

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

    }; // class Ensemble

  } // namespace assemble
} // namespace bcl

#endif // BCL_ASSEMBLE_ENSEMBLE_H_

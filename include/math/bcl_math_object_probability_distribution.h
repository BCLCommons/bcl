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

#ifndef BCL_MATH_OBJECT_PROBABILITY_DISTRIBUTION_H_
#define BCL_MATH_OBJECT_PROBABILITY_DISTRIBUTION_H_

// include the namespace header
#include "bcl_math.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"
#include "util/bcl_util_implementation.h"
#include "util/bcl_util_serializable_interface.h"
#include "util/bcl_util_sh_ptr_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ObjectProbabilityDistribution
    //! @brief Class for storing objects and corresponding probabilities in one architecture
    //!
    //! @see @link example_math_object_probability_distribution.cpp @endlink
    //! @author woetzen, karakam, fischea
    //! @date Nov 20, 2006
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    template< typename t_ObjectType>
    class ObjectProbabilityDistribution :
      public util::SerializableInterface
    {

    ///////////
    // types //
    ///////////

    public:

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //!
      //! @class Assignment
      //! @brief stores an object and its corresponding probability
      //!
      //! @remarks example unnecessary
      //! @author fischea
      //! @date Oct 07, 2016
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      class Assignment :
        public storage::Pair< double, util::Implementation< t_ObjectType> >,
        public util::SerializableInterface
      {

      private:

        typedef storage::Pair< double, util::Implementation< t_ObjectType> > Base;

      public:

        //! @brief default constructor
        Assignment() :
          Base()
        {
        }

        //! @brief construct Pair from two values
        //! @param FIRST the object which will be m_First
        //! @param SECOND the object which will be m_Second
        Assignment
        (
          const double &FIRST,
          const t_ObjectType &SECOND
        ) :
          Base( FIRST, SECOND)
        {
        }

        //! @brief construct Pair from std::pair
        //! @param PAIR which will be the Pair
        Assignment( const Base &PAIR) :
          Base( PAIR)
        {
        }

        //! copy constructor
        Assignment( const Assignment &PAIR) :
          Base( PAIR)
        {
        }

        //! @brief clone function
        //! @return pointer to a new Assignment
        Assignment *Clone() const
        {
          return new Assignment( *this);
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

        //! @brief alias
        //! @return an alias for the class
        const std::string &GetAlias() const
        {
          return this->Second().GetAlias();
        }

        //! @brief alias
        //! @return an alias for the class
        const std::string &GetScheme() const
        {
          static const std::string s_empty;
          return s_empty;
        }

        //! @brief return parameters for member data that are set up from the labels
        //! @return parameters for member data that are set up from the labels
        io::Serializer GetSerializer() const
        {
          io::Serializer parameters;
          parameters.SetClassDescription( "object and corresponding probability");
          parameters.AddInitializer
          (
            "probability",
            "probability of the object",
            io::Serialization::GetAgent( &this->First())
          );
          parameters.AddInitializer
          (
            "",
            "the object of interest",
            io::Serialization::GetAgent( &this->Second())
          );

          return parameters;
        }

      };

    ////////
    //data//
    ////////

    private:

      //! List to store probabilities and associated objects to be returned
      storage::Vector< Assignment> m_ProbabilityObjectList;

    //////////
    // data //
    //////////

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      ObjectProbabilityDistribution()
      {
      }

      //! @brief constructor
      ObjectProbabilityDistribution
      (
        const storage::Vector< Assignment> &PROBABILITY_OBJECT_LIST
      ) :
        m_ProbabilityObjectList( PROBABILITY_OBJECT_LIST)
      {
      }

      //! @brief constructor from an object probability distribution of another type
      //! @tparam t_OtherData data type in parameter. must be able to create a t_ObjectType from t_OtherData
      //! @param PROBABILITY_OBJECT probability object that will be used to create this
      template< typename t_OtherData>
      ObjectProbabilityDistribution
      (
        const ObjectProbabilityDistribution< t_OtherData> &PROBABILITY_OBJECT
      ) :
        m_ProbabilityObjectList()
      {
        // iterate through the PROBABILITY_OBJECT
        for
        (
          auto itr( PROBABILITY_OBJECT.GetData().Begin()), itr_end( PROBABILITY_OBJECT.GetData().End());
          itr != itr_end;
          ++itr
        )
        {
          // make a t_ObjectType from a t_OtherData
          const t_ObjectType &this_datatype( *( *itr).Second());

          // add this datatype to m_ProbabilityObjectList
          m_ProbabilityObjectList.PushBack( ( *itr).First(), this_datatype);
        }
      }

      //! @brief virtual copy constructor
      virtual ObjectProbabilityDistribution< t_ObjectType> *Clone() const
      {
        return new ObjectProbabilityDistribution< t_ObjectType>( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      virtual const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief gives the vector of probabilities and objects
      //! @return the vector of probabilities and objects which is m_ProbabilityObjectList
      const storage::Vector< Assignment> &GetData() const
      {
        return m_ProbabilityObjectList;
      }

      //! @brief get a short name for this class
      //! @return a short name for this class
      const std::string &GetAlias() const
      {
        static const std::string s_name( "ObjectProbabilityDistribution");
        return s_name;
      }

      //! @brief returns parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const
      {
         io::Serializer serializer;
         serializer.SetClassDescription( "Stores objects with corresponding probabilities.");
         serializer.AddInitializer
         (
           "probabilities",
           "objects with corresponding probabilities",
           io::Serialization::GetAgent( &m_ProbabilityObjectList)
         );

         return serializer;
      }

    ///////////////
    // operators //
    ///////////////

      //! @brief push back the given pair of probability of object and pair
      //! @param PROBABILITY_OBJECT_PAIR pair of object and its according probability
      void PushBack( const Assignment &PROBABILITY_OBJECT_PAIR)
      {
        m_ProbabilityObjectList.PushBack( PROBABILITY_OBJECT_PAIR);
      }
      //! @brief push back the given pair of probability of object and pair
      //! @param PROBABILITY, OBJECT pair of object and its according probability
      void PushBack( const double &PROBABILITY, const t_ObjectType &OBJECT)
      {
        m_ProbabilityObjectList.PushBack( Assignment( PROBABILITY, OBJECT));
      }

      //! @brief Return a random object according to the supplied probabilities
      //! @return a randomly drawn object from the vector
      const t_ObjectType &DetermineRandomCase() const
      {
        // determining the random value
        double temp( random::GetGlobalRandom().Random< double>( CalculateSum()));

        for
        (
          auto prob_itr( m_ProbabilityObjectList.Begin()), prob_itr_end( m_ProbabilityObjectList.End());
          prob_itr != prob_itr_end;
          ++prob_itr
        )
        {
          temp -= ( prob_itr)->First();
          if( temp <= 0)
          {
            return *( prob_itr)->Second();
          }
        }

        BCL_Exit( "The probabilities are not set correctly or distribution is empty!", -1);
        return *( *m_ProbabilityObjectList.End()).Second();
      }

      //! @brief Calculates the sum of probabilities
      //! @return the sum of all probabilities in the vector
      double CalculateSum() const
      {
        // instantiate sum
        double sum( 0);

        // iterate over all probabilities and sum them up
        for
        (
          auto prob_itr( m_ProbabilityObjectList.Begin()), prob_itr_end( m_ProbabilityObjectList.End());
          prob_itr != prob_itr_end;
          ++prob_itr
        )
        {
          sum += ( prob_itr)->First();
        }
        // return sum
        return sum;
      }

      //! @brief normalizes the variable to sum
      //! @param SUM requested sum of the distribution weights after normalization
      void SetToSum( const double SUM = 1.0)
      {
        // calculate and store ratio
        const double ratio( SUM / CalculateSum());

        // iterate over all pairs of weights and objects
        for
        (
          auto prob_itr( m_ProbabilityObjectList.Begin()), prob_itr_end( m_ProbabilityObjectList.End());
          prob_itr != prob_itr_end;
          ++prob_itr
        )
        {
          // update the weight
          ( *prob_itr)->First() *= ratio;
        }
      }

    //////////////////////
    // input and output //
    //////////////////////

    }; // template class ObjectProbabilityDistribution

    // instantiate s_Instance
    template< typename t_ObjectType>
    const util::SiPtr< const util::ObjectInterface> ObjectProbabilityDistribution< t_ObjectType>::s_Instance
    (
      util::Enumerated< t_ObjectType>::AddInstance( new ObjectProbabilityDistribution< t_ObjectType>())
    );

  } // namespace math
} // namespace bcl

#endif //BCL_MATH_OBJECT_PROBABILITY_DISTRIBUTION_H_

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

#ifndef BCL_ASSEMBLE_TOPOLOGY_DISTANCE_H_
#define BCL_ASSEMBLE_TOPOLOGY_DISTANCE_H_

// include the namespace header
#include "bcl_assemble.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "util/bcl_util_object_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class TopologyDistance
    //! @brief class is a measure to compare the topologies of two models for the same protein.
    //! @details This measure takes in two ProteinModels. assuming one is a native/template model and they are already
    //! superimposed on each other, and then evaluates the superimposition by looking at the whether the same SSEs
    //! exists, if they are in similarly same locations and with the same orientations.
    //!
    //! @see @link example_assemble_topology_distance.cpp @endlink
    //! @author karakam
    //! @date Oct 6, 2009
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API TopologyDistance :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      //! Minimum distance between the centers of two corresponding SSEs in separate models to be considered significantly different
      double m_DistanceCutoff;

      //! Minimum angle change between the orientations of two corresponding SSEs in separate models to be considered significantly different
      double m_AngleCutoff;

    public:

    //////////
    // data //
    //////////

      //! @brief returns default distance cutoff
      //! @return default distance cutoff
      static double GetDefaultDistanceCutoff();

      //! @brief returns default angle cutoff
      //! @return default angle cutoff
      static double GetDefaultAngleCutoff();

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      TopologyDistance();

      //! @brief constructor from a distance cutoff and angle cutoff
      //! @param DISTANCE_CUTOFF distance cutoff
      //! @param ANGLE_CUTOFF angle cutoff
      TopologyDistance( const double DISTANCE_CUTOFF, const double ANGLE_CUTOFF);

      //! @brief Clone function
      //! @return pointer to new TopologyDistance
      TopologyDistance *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief returns distance cutoff
      //! @return distance cutoff
      double GetDistanceCutoff() const;

      //! @brief returns angle cutoff
      //! @return angle cutoff
      double GetAngleCutoff() const;

    ////////////////
    // operations //
    ////////////////

    ///////////////
    // operators //
    ///////////////

      //! @brief calculates the distance between two given protein models
      //! @param MODEL ProteinModel for which distance to given TEMPLATE_MODEL will be calculated
      //! @param TEMPLATE_MODEL ProteinModel that will be used as a template
      double operator()( const ProteinModel &MODEL, const ProteinModel &TEMPLATE_MODEL) const;

      //! @brief calculates the distance between two given chains
      //! @param CHAIN Chain for which distance to given TEMPLATE_CHAIN will be calculated
      //! @param TEMPLATE_CHAIN Chain that will be used as a template
      double operator()( const Chain &CHAIN, const Chain &TEMPLATE_CHAIN) const;

      //! @brief calculates the distance between two given SSEs
      //! @param SS_ELEMENT for which distance to given TEMPLATE_SS_ELEMENT will be calculated
      //! @param TEMPLATE_SS_ELEMENT SSE that will be used as a template
      double operator()( const SSE &SS_ELEMENT, const SSE &TEMPLATE_SS_ELEMENT) const;

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

    //////////////////////
    // helper functions //
    //////////////////////

    private:

    }; // class TopologyDistance

  } // namespace assemble
} // namespace bcl

#endif // BCL_ASSEMBLE_TOPOLOGY_DISTANCE_H_ 

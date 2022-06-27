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

#ifndef BCL_CHEMISTRY_CONFORMATION_INTERFACE_H_
#define BCL_CHEMISTRY_CONFORMATION_INTERFACE_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically
#include "storage/bcl_storage.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from the bcl
#include "bcl_chemistry_atom_complete.h"
#include "bcl_chemistry_atom_conformational_shared.h"
#include "bcl_chemistry_configuration_interface.h"
#include "bcl_chemistry_constitution_interface.h"
#include "bcl_chemistry_has_properties_interface.h"
#include "coord/bcl_coord_movable_interface.h"
#include "descriptor/bcl_descriptor_has_cache.h"
#include "descriptor/bcl_descriptor_sequence_interface.h"
#include "graph/bcl_graph_undirected_edge.h"
#include "iterate/bcl_iterate_generic.h"
#include "linal/bcl_linal_vector_3d.h"
#include "sdf/bcl_sdf_atom_info.h"
#include "signal/bcl_signal_signal.h"
#include "storage/bcl_storage_triplet.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ConformationInterface
    //! @brief Interface class for molecule/fragment conformation
    //! @details Interface class for conformation of molecule or fragment. Handles molecular/fragment conformation data.
    //!
    //! @remarks example unnecessary
    //! @author kothiwsk, mendenjl, brownbp1
    //! @date Dec 08, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ConformationInterface :
      public descriptor::SequenceInterface< AtomConformationalInterface>,
      public HasPropertiesInterface< coord::MovableInterface>
    {
    //////////
    // data //
    //////////

      //! Change signal for conformation
      mutable signal::Signal1< const ConformationInterface &> m_ChangeSignal;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief desctruction for conformation
      virtual ~ConformationInterface()
      {
        // emit the change signal
        m_ChangeSignal.Emit( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns an iterator of atom conformations
      //! @return a generic iterator for atom conformations
      virtual iterate::Generic< const AtomConformationalInterface> GetAtomsIterator() const = 0;

      //! @brief returns an iterator of atom conformations
      //! @return a generic iterator for atom conformations
      virtual iterate::Generic< const AtomConformationalInterface> GetIterator() const
      {
        return GetAtomsIterator();
      }

      //! @brief access to the change signal handler
      //! @return const ref to the SignalHandler that emits on changes
      signal::Signal1< const ConformationInterface &> &GetChangeSignal() const
      {
        return m_ChangeSignal;
      }

      //! @brief return the number of atoms
      //! @return the number of atoms
      virtual const size_t &GetNumberAtoms() const = 0;

      //! @brief return the number of bonds
      //! @return the number of bonds
      virtual const size_t &GetNumberBonds() const = 0;

      //! @brief return a number of hydrogens in the conformation
      //! @return the number of hydrogens in the conformation
      size_t GetNumberHydrogens() const;

      //! @brief return the number of bonds with a particular property (e.g. is in ring, is aromatic, is isometric)
      //! @param DATA the data to retrieve for each bond
      //! @param VALUE the value to count, e.g. if ConfigurationalBondType->GetData( DATA) == VALUE, ++return
      //! @return the number of bonds with the property of interest
      size_t CountNonValenceBondsWithProperty
      (
        const ConfigurationalBondTypeData::Data &DATA,
        const size_t &VALUE
      ) const;

      //! @brief return the adjacency list
      //! @param BOND_SCHEME how to represent bond data as a size_t
      //! @return the adjacency list
      virtual storage::Vector< graph::UndirectedEdge< size_t> >
        GetAdjacencyList( const ConfigurationalBondTypeData::Data &BOND_SCHEME) const = 0;

      //! @brief get bonds of molecule
      //! @return all the connectivities of a molecule
      virtual storage::Vector< sdf::BondInfo> GetBondInfo() const = 0;

      //! @brief get info on all atoms in the molecule
      //! @return all the atom info about the molecule
      virtual storage::Vector< sdf::AtomInfo> GetAtomInfo() const = 0;

      //! @brief get the index of a particular atom conformational interface
      //! @param ATOM the atom of interest
      //! @return the index of a particular atom conformational interface (undefined if atom is not in this molecule)
      virtual size_t GetAtomIndex( const AtomConformationalInterface &ATOM) const = 0;

      //! @brief get the number of valences for the entire molecule or fragment
      //! @return the number of valences on the entire molecule or fragment
      size_t GetNumberValences() const;

      //! @return a vector with the atom types in the molecule
      storage::Vector< AtomType> GetAtomTypesVector() const;

      //! @return a string with each atom type in the molecule
      std::string GetAtomTypesString() const;

      //! @return a string with each element type in the molecule
      std::string GetElementTypesString() const;

      //! @return a string with each bond type in the molecule
      std::string GetBondTypesString() const;

      //! @brief return a string containing the chirality of the atoms of this molecule as a string
      std::string GetChiralityString() const;

      //! @brief return a string containing the double bond isometry of the atoms of this molecule as a string
      std::string GetDoubleBondIsometryString() const;

      //! @brief returns a string with the sum formula
      //! @return a string with the sum formula
      std::string GetSumFormula() const;

      //! @brief assignment operator
      virtual ConformationInterface &operator=( const ConformationInterface &C);

      //! @brief checks for 'bad' geometry, i.e., z-coord./all coords. == const.
      //! @return flag whether geometry is 'bad'
      storage::Vector< size_t> GetAtomsWithBadGeometry() const;

      //! @brief checks for 'bad' geometry, i.e., z-coord./all coords. == const.
      //! @return flag whether geometry is 'bad'
      bool HasBadGeometry() const;

      //! @brief checks for non-gasteiger atom types
      //! @return flag whether any non-gasteiger types are present
      bool HasNonGasteigerAtomTypes() const;

      //! @brief checks for chiral centers
      //! @return flag whether any chiral centers are present
      bool HasChiralCenters() const;

      //! @brief checks for double bond isometry
      //! @return double bond isometry
      bool HasIsometry() const;

      //! @brief checks any arbitrary generic iterator over AtomConformationInterfcaes for non-gasteiger atom types
      //! @param ATOMS a generic iterator over the atoms
      //! @return flag whether any non-gasteiger types are present
      static bool HasNonGasteigerAtomTypes( iterate::Generic< const AtomConformationalInterface> ATOMS);

      //! @brief get the center of the molecule, as defined by the positions
      linal::Vector3D GetCenter() const;

      //! @brief return a SiPtrVector to linal::Vector3D containing atomic coordinates
      util::SiPtrVector< const linal::Vector3D> GetAtomCoordinates() const;

      //! @brief return a SiPtrVector to linal::Vector3D containing heavy atom coordinates
      util::SiPtrVector< const linal::Vector3D> GetHeavyAtomCoordinates() const;

      //! @brief calculate the bond angles
      //! @return a map from atom conformation interface in ascending order of atom position
      storage::Vector< double> GetBondAngles() const;

      //! @brief calculate all the dihedral angles
      //! @return a vector containing dihedral angles, in ascending order of atom positions
      storage::Vector< double> GetDihedralAngles() const;

      //! @brief calculate all the dihedral angles
      //! @return a vector containing dihedral angles, in ascending order of atom position
      storage::Map< storage::VectorND< 7, size_t>, storage::Vector< double> > GetDihedralAnglesByType() const;

      //! @brief calculate all the dihedral angles
      //! @return a vector containing dihedral angles, in ascending order of atom positions
      storage::Map< storage::VectorND< 5, size_t>, storage::Vector< double> > GetBondAnglesByType() const;

      //! @brief calculate all the bond lengths
      //! @return a vector containing bond lengths
      storage::Map< storage::VectorND< 3, size_t>, storage::Vector< double> > GetBondLengthsByType() const;

      //! @brief get name of small molecule
      virtual const std::string &GetName() const = 0;

      //! @brief set name of small molecule
      virtual void SetName( const std::string &NAME) = 0;

      //! @brief get the conformational properties merged with the constitutional and configurational properties
      //! if a properties are in multiple layers, the conformation's value will be used
      virtual SmallMoleculeMiscProperties GetMergedProperties() const = 0;

      //! @brief check to see if amide bonds are planer
      //! @param TOLERANCE tolerance (in degrees) for non-planarity
      //! @return true if bonds are planer else return false
      bool AreAmideBondsPlaner( const double &TOLERANCE) const;

      //! @brief compute deviations from planarity for all amide bonds
      //! @return a vector containing the indices of the N and C bonded amide atoms and
      //! the corresponding deviation from nonplanarity
      //! @return a vector of triplets where the vector is indexed by the amide bond and the triplet contains
      //! atom indices corresponding to the lower and higher amide bond atom indices, respectively, and
      //! the magnitude of the deviation
      storage::Vector< storage::Triplet< size_t, size_t, double> > GetAmideBondNonPlanarity() const;

      //! @brief Sum the deviation of amide bond planarity over the molecule
      //! @param AMIDE_DEVIATIONS from GetAmideBondNonPlanarity()
      //! @param TOLERANCES (in degrees) acceptable deformations of the amide
      //! @param PENALTIES (in BCL::Conf units) for amide deviations exceeding the tolerance;
      //! must be same size as TOLERANCES
      //! @return the per-amide bond penalty
      storage::Vector< double> GetPerAmideBondNonPlanarityPenalty
      (
        storage::Vector< storage::Triplet< size_t, size_t, double> > &AMIDE_DEVIATIONS,
        storage::Vector< double> &TOLERANCES,
        storage::Vector< double> &PENALTIES
      ) const;

      //! @brief Sum the deviation of amide bond planarity over the molecule
      //! @param TOLERANCE tolerance (in degrees) for non-planarity, amide bonds with deviation smaller than tolerance
      //!        will contribute 0 to the sum. 25 degrees is recommended, because highly-strained amide bonds (e.g. in
      //!        tetramethylurea) can deviate up to 22 degrees from planar due to steric stress.
      //! @return deviation of amide bond planarity summed over the molecule
      double GetTotalAmideBondNonPlanarity( const double &TOLERANCE) const;

      //! @brief check to see if aromatic rings are planer
      //! @return true if aromatic rings are planar, if not return false
      bool AreAromaticRingsPlaner() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief move the molecule
      //! @param TRANSLATION amount to change x,y,z coordinates
      void Translate( const linal::Vector3D &TRANSLATION);

      //! @brief transform the object by a given TransformationMatrix3D
      //! @param TRANSFORMATION_MATRIX_3D TransformationMatrix3D to be applied
      void Transform( const math::TransformationMatrix3D &TRANSFORMATION_MATRIX_3D);

      //! @brief rotate the object by a given RotationMatrix3D
      //! @param ROTATION_MATRIX_3D RotationMatrix3D to be applied
      void Rotate( const math::RotationMatrix3D &ROTATION_MATRIX_3D);

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @return output stream which was written to
      std::ostream &WriteMDL( std::ostream &OSTREAM) const;

    private:

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief get changable atoms conformational interface
      //! @return iterator to changable atoms conformational interface, which allows this class to call SetPosition
      virtual iterate::Generic< AtomConformationalInterface> GetAtomsIteratorNonConst() = 0;

    }; // class ConformationInterface

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_CONFORMATION_INTERFACE_H_

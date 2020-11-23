/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.
    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.
    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "circularDeformation.H"
#include "pointPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "Time.H"
#include "polyMesh.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

circularDeformation::
circularDeformation
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchField<vector>(p, iF),
    b_(0.0),
    epsilon_(0.0),
    omg_(0.0),
    kmaxscale_(0.0),
    p0_(p.localPoints())
{}


circularDeformation::
circularDeformation
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const dictionary& dict
)
:
    fixedValuePointPatchField<vector>(p, iF, dict),
    b_(readScalar(dict.lookup("b"))),
    epsilon_(readScalar(dict.lookup("eps"))),
    omg_(readScalar(dict.lookup("omg"))),
    kmaxscale_(readScalar(dict.lookup("k")))
{
    if (!dict.found("value"))
    {
        updateCoeffs();
    }

    if (dict.found("p0"))
    {
        p0_ = vectorField("p0", dict , p.size());
    }
    else
    {
        p0_ = p.localPoints();
    }
}


circularDeformation::
circularDeformation
(
    const circularDeformation& ptf,
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
    fixedValuePointPatchField<vector>(ptf, p, iF, mapper),
    b_(ptf.b_),
    epsilon_(ptf.epsilon_),
    omg_(ptf.omg_),
    kmaxscale_(ptf.kmaxscale_),
    p0_(ptf.p0_, mapper)
{}


circularDeformation::
circularDeformation
(
    const circularDeformation& ptf,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchField<vector>(ptf, iF),
    b_(ptf.b_),
    epsilon_(ptf.epsilon_),
    omg_(ptf.omg_),
    kmaxscale_(ptf.kmaxscale_),
    p0_(ptf.p0_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void circularDeformation::autoMap
(
    const pointPatchFieldMapper& m
)
{
    fixedValuePointPatchField<vector>::autoMap(m);

    p0_.autoMap(m);
}


void circularDeformation::rmap
(
    const pointPatchField<vector>& ptf,
    const labelList& addr
)
{
    const circularDeformation& oVptf =
        refCast<const circularDeformation>(ptf);

    fixedValuePointPatchField<vector>::rmap(oVptf, addr);

    p0_.rmap(oVptf.p0_, addr);
}


void circularDeformation::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const polyMesh& mesh = this->internalField().mesh()();
    const Time& t = mesh.time();
    const pointPatch& p = this->patch();
    vectorField sd = p0_;
    

    scalar kmax = epsilon_*kmaxscale_/b_;
    scalar A = epsilon_*b_;
 

    scalar a1 = A*sin(omg_*t.value());
    scalar k1 = kmax*sin(omg_*t.value());
    scalar rad = 1/(VSMALL+k1);
    scalar theta = b_/2/rad;

    
    forAll(p0_,iter)
    {
        sd[iter] = vector(rad*sin(theta*p0_[iter][0]*2/b_)-p0_[iter][0], (rad*cos(theta*p0_[iter][0]*2/b_)-rad+a1), 0);
    };

    Field<vector>::operator=
    (
        (p0_+sd-p.localPoints())
       /t.deltaTValue()
    );

    fixedValuePointPatchField<vector>::updateCoeffs();
}


void circularDeformation::write(Ostream& os) const
{
    const pointPatch& p = this->patch();
    pointPatchField<vector>::write(os);
    os.writeKeyword("b")
        << b_ << token::END_STATEMENT << nl;
    os.writeKeyword("eps")
        << epsilon_ << token::END_STATEMENT << nl;
        os.writeKeyword("omg")
        << omg_ << token::END_STATEMENT << nl;
        os.writeKeyword("k")
        << kmaxscale_ << token::END_STATEMENT << nl;
    p0_.writeEntry("p0", os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePointPatchTypeField
(
    pointPatchVectorField,
    circularDeformation
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

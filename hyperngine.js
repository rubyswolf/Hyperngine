const tau = Math.PI*2 //Define the superior circle constant

class Complex {
    constructor(real, imag) {
        this.real = real;
        this.imag = imag;
    }

    add(other) {
        return new Complex(this.real + other.real, this.imag + other.imag);
    }

    subtract(other) {
        return new Complex(this.real - other.real, this.imag - other.imag);
    }

    get conjugate() {
        return new Complex(this.real, -this.imag);
    }

    get squareMag() {
        return this.real ** 2 + this.imag ** 2;
    }

    multiply(other) {
        return new Complex(
            this.real * other.real - this.imag * other.imag,
            this.real * other.imag + this.imag * other.real
        );
    }

    dot(other) { //Dot product of two complex numbers
        return this.real * other.real + this.imag * other.imag
    }

    divide(other) {
        const denom = other.real ** 2 + other.imag ** 2;
        return new Complex(
            (this.real * other.real + this.imag * other.imag) / denom,
            (this.imag * other.real - this.real * other.imag) / denom
        );
    }

    get asHyperbolicVector() {
        return new HyperbolicVector(this.real,this.imag)
    }
}

class HyperbolicVector { //A vector in the poincare disk model
    constructor(x, y) {
        this.x = x;
        this.y = y;
    }

    static fromPolar(r, theta) { //Create a Hyperbolic Vector from polar form
        const poincareR = Math.tanh(r/2); //Convert from true magnitude to poincare disk magnitude
        return new HyperbolicVector(
            poincareR * Math.cos(theta),
            poincareR * Math.sin(theta)
        );
    }

    add(other) { //Adds `other` to `this` tip to tail. Note that in hyperbolic geometry u+v is not the same as v+u
        let u = this.asComplex
        let v = other.asComplex
        
        let norm1Sq=new Complex(u.squareMag,0)
        let norm2Sq=new Complex(v.squareMag,0)
        let uv2=new Complex(2*u.dot(v),0)
        return (((new Complex(1,0).add(uv2).add(norm2Sq)).multiply(u)).add((new Complex(1,0).subtract(norm1Sq)).multiply(v))).divide(new Complex(1,0).add(uv2).add(norm1Sq.multiply(norm2Sq))).asHyperbolicVector
    }

    get asComplex() {
        return new Complex(this.x,this.y)
    }

    get asPolar() {
        return new PolarHyperbolicVector(this.magnitude,this.angle)
    }

    to(other) { //Finds the vector from `this` to `other` such that this+result=other
        let v = this.asComplex
        let x = other.asComplex
        return (x.subtract(v)).divide(new Complex(1,0).subtract(v.conjugate.multiply(x))).asHyperbolicVector
    }

    scale(factor) {           
        return HyperbolicVector.fromPolar(this.magnitude * factor, this.angle);
    }

    distance(other) { //Finds the true distance from `this` to `other`
        return this.to(other).magnitude
    }

    angleBetween(other1, other2) //Finds the unsigned angle at `this` between the two geodesics from `this` to `other1` and from `this` to `other2`
    {
        return (this.to(other1).angle-this.to(other2).angle+tau)%tau
    }

    rotate(rotation) {
        return HyperbolicVector.fromPolar(this.magnitude,this.angle+rotation)
    }

    project(projection) //Projects the vector to a different projection from a different origin with a different rotation
    {
        function projectTanh(x,t) { //Just y=tanh(x) when t=0 and just y=x when t=1, used to interpolate between stereographic and orthographic
            if (t==1)
            {
                return x
            }
            let t2=1/(1-t*t)
            return Math.tanh(x/t2)*t2
        }

        let shifted=projection.position.to(this.scale(projection.scale)) //Shift the incoming vector to the new origin
        const r = projectTanh(((shifted.magnitude)/projection.factor),projection.orthographic); //The new radius from the centre as defined by the projection
        const theta = shifted.angle+projection.rotation //Adds an extra rotation
        return new ProjectedVector( //Convert to a projected vector instance
            r * Math.cos(theta),
            r * Math.sin(theta)
        );
    }

    subject(subRadius) { //Projects the vector into the poincare model a second time
        let poincareR = Math.tanh((this.diskMagnitude*subRadius)/2);
        let theta = this.angle 
        return new HyperbolicVector(
            poincareR * Math.cos(theta),
            poincareR * Math.sin(theta)
        );
    }

    reflect(geo) { //Reflects a point across a geodesic
        let cent = geo.centre
        let r = geo.radius
        let dist = (this.x-cent.x)**2+(this.y-cent.y)**2
        return new HyperbolicVector(cent.x+((r*r)/dist)*(this.x-cent.x),cent.y+((r*r)/dist)*(this.y-cent.y))
    }
    
    get diskMagnitude() { //The magnitude of the point in the disk
        return Math.sqrt(this.x ** 2 + this.y ** 2);
    }

    get magnitude() { //True true hyperbolic distance from the origin
        return 2*Math.atanh(this.diskMagnitude);
    }

    get angle() {
        return Math.atan2(this.y, this.x);
    }
}

class PolarHyperbolicVector { //A hyperbolic vector in polar form
    constructor (magnitude, angle) {
        this.magnitude = magnitude
        this.angle = angle
    }

    rotate(rotation) {
        return new PolarHyperbolicVector(this.magnitude, this.angle+rotation)
    }

    scale(factor) {           
        return new PolarHyperbolicVector(this.magnitude*factor, this.angle)
    }

    get asHyperbolicVector() {
        return HyperbolicVector.fromPolar(this.magnitude,this.angle)
    }
}

class Geodesic { //Defines a geodesic (The shortest path between two points which you can think of as essentially just a straight line)
    constructor (a,b) {
        this.a=a //First point
        this.b=b //Second point
    }

    //In the poincare model, geodesics are represended as arcs of a circle that meet the boundry of the disk at a right angle

    get centre() { //Centre of the geodisc arc
        let x1=this.a.x
        let y1=this.a.y
        let x2=this.b.x
        let y2=this.b.y
        let denom=(x1*x1+y1*y1)
        let x3=x1/denom
        let y3=y1/denom
        let m1=(x1-x2)/(y2-y1)
        let c1=(y1+y2+((x2-x1)/(y2-y1))*(x1+x2))/2
        let m2 = (y3-y1)==0 ? 0 : (x1-x3)/(y3-y1)
        let c2=(y1+y3-m2*(x1+x3))/2
        let x4=y3-y1==0?(x1+x3)/2:(c2-c1)/(m1-m2)
        let y4=m1*x4+c1
        return new HyperbolicVector(x4,y4)
    }
    
    get radius() //Radius of the geodesic arc
    {
        let cent = this.centre
        return Math.sqrt((cent.x-this.a.x)**2+(cent.y-this.a.y)**2)
    }

    get midpoint() { //Midpoint of the geodesic
        return this.along(0.5)
    }

    along(t) { //A point along the geodesic from 0 being a and 1 being b
        return this.a.add(this.a.to(this.b).scale(t))
    }
}

class ProjectedVector { //A vector that has been projected from poincare
    constructor(x,y) {
        this.x=x;
        this.y=y;
    }
    
    get diskMagnitude() {
        return Math.sqrt(this.x ** 2 + this.y ** 2);
    }
    
    get magnitude() {
        return 2*Math.atanh(this.diskMagnitude);
    }

    get angle() {
        return Math.atan2(this.y, this.x);
    }

    unproject(projection) { //The same as project for HyperbolicVector but everything happens in reverse to restore the original
        function unprojectTanh(x,t) {
            if (t==1)
            {
                return x
            }
            let t2=1/(1-t*t)
            return Math.atanh(x/t2)*t2
        }

        const r = Math.tanh((unprojectTanh(this.diskMagnitude,projection.orthographic)*projection.factor)/2);
        const theta = this.angle-projection.rotation
        return projection.position.add(new HyperbolicVector(
            r * Math.cos(theta),
            r * Math.sin(theta)
        )).scale(1/projection.scale);
    }
}

class Projection {
    constructor(position,rotation,scale,factor,orthographic) {
        //Projection controls the distance from the disk to the hyperboloid in stereographic mode.
        //projection=1 gives the Beltrami-Klein model
        //projection=2 gives the Poincare model
        //Anything inbetween interpolates betwen them
        //In orthographic mode it basically just controls the zoom factor

        //Rotation adds a rotational offset

        //Orthographic controls the projection type from stereographic to orthographic
        //orthographic=0 gives stereographic projections such as Klein and Poincare
        //orthographic=1 gives a special kind of projection where distance from the origin when projection=1 is the true distance to that point
        //This not exactly the same as a true orthographic projection by raycasting from the hyperboloid to the plane

        //Position controls where the new origin should be

        this.position=position
        this.rotation=rotation
        this.scale=scale
        this.factor=factor
        this.orthographic=orthographic
    }
}

class ShvgRenderer { //A renderer that can render a hyperbolic scene into an SVG
    constructor(projection) {
        this.projection=projection
        this.scale=1
        this.calculateProperties()
    }

    get size() {
        return this.baseSize*this.scale
    }
    
    calculateProperties() {
        this.baseSize = Math.sqrt(window.innerWidth**2+window.innerHeight**2)/2;

        this.centreX = window.innerWidth / 2;
        this.centreY = window.innerHeight / 2;
        
        this.scale = Math.max(Math.min(this.scale,1),Math.min(window.innerWidth,window.innerHeight)/Math.sqrt(window.innerWidth**2+window.innerHeight**2))
    }

    toSVGCoords(vector) {
        vector = vector.project(this.projection)
        return {
            x: this.centreX + vector.x * this.size,
            y: this.centreY - vector.y * this.size
        };
    }

    generateBaseSVG() {
        return `<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${window.innerWidth} ${window.innerHeight}">
            <circle cx="${this.centreX}" cy="${this.centreY}" r="${this.size}" 
                    fill="#0B2126" stroke="none" stroke-width="1"/>`;
    }

    addPoint(vector, colour = "blue", radius = 7) {
        const {x, y} = this.toSVGCoords(vector);
        return `<circle cx="${x}" cy="${y}" r="${radius}" fill="${colour}"/>`;
    }

    geoPath(geo,samples,continuing = false) {
        let points = [];
        
        // Generate sample points along the geodesic
        for (let i = 0; i < samples; i++) {
            let t = i / (samples - 1);
            let point = geo.along(t);
            points.push(this.toSVGCoords(point));
        }

        // Start path from the first point
        let path = continuing?``:`M ${points[0].x} ${points[0].y} `

        // Convert points along the geodesic to Bézier curves
        for (let i = 0; i < points.length - 1; i++) {
            let p0 = points[Math.max(i - 1, 0)];
            let p1 = points[i];
            let p2 = points[i + 1];
            let p3 = points[Math.min(i + 2, points.length - 1)];

            // Calculate control points for cubic Bézier from Catmull-Rom spline
            let controlPoint1 = {
                x: p1.x + (p2.x - p0.x) / 6,
                y: p1.y + (p2.y - p0.y) / 6
            };
            let controlPoint2 = {
                x: p2.x - (p3.x - p1.x) / 6,
                y: p2.y - (p3.y - p1.y) / 6
            };

            // Use cubic Bézier curve to approximate spline segment
            path += `C ${controlPoint1.x} ${controlPoint1.y}, ${controlPoint2.x} ${controlPoint2.y}, ${p2.x} ${p2.y} `;
        }
        return path
    }

    addGeodesic(geo,stroke,samples) {
        return `<path d="${this.geoPath(geo,samples)}" stroke="${stroke.colour}" stroke-width="${stroke.width}" fill="none"/>`;
    }

    addPolygon(points,samples,fill,stroke) {
        let path = points.map(((e,i)=>{return this.geoPath(new Geodesic(e,points[(i+1+points.length)%points.length]),samples,i!=0)})).join("")+"Z"
        return stroke.enabled?`<path d="${path}" stroke="${stroke.colour}" stroke-width="${stroke.width}" fill="${fill}"/>`:`<path d="${path}" fill="${fill}"/>`;
    }

    addBlob(points,fill,stroke) {
        points = points.map(point=>{return this.toSVGCoords(point)})
        if (points.length < 2) return ''; // Ensure there are enough points to form a shape

        // Calculate the initial midpoint between the first and last points
        const firstMidpoint = {
            x: (points[0].x + points[points.length - 1].x) / 2,
            y: (points[0].y + points[points.length - 1].y) / 2
        };
        
        // Start the path from the first midpoint
        let path = `M ${firstMidpoint.x},${firstMidpoint.y}`;
        
        // Loop through each point and create quadratic Bézier segments to the next midpoint
        for (let j = 0; j < points.length; j++) {
            const nextIndex = (j + 1) % points.length; // Wrap around to form a closed shape
            const midpoint = {
                x: (points[j].x + points[nextIndex].x) / 2,
                y: (points[j].y + points[nextIndex].y) / 2
            };
            path += ` Q ${points[j].x},${points[j].y} ${midpoint.x},${midpoint.y}`;
        }
    
        // Return the complete SVG path element
        return stroke.enabled?`<path d="${path}" fill="${fill}" stroke="${stroke.colour}" stroke-width="${stroke.width}"/>`:`<path d="${path}" fill="${fill}""/>`;
    }

    addCircle(centre,radius,samples,fill,stroke=Stroke.none) {
        return this.addBlob(Array.from({ length: samples }, (v, i) => {return centre.add(HyperbolicVector.fromPolar(radius,(i/samples)*tau))}),fill,stroke)
    }

    render(elements) {
        let svg = this.generateBaseSVG();
        svg += elements.join('\n');
        svg += '</svg>';
        return svg;
    }
}

class Stroke {
    constructor(colour,width,fixed,enabled=true) {
        this.enabled=enabled
        if (enabled)
        {
            this.colour=colour
            this.width=width
            this.fixed=fixed
        }
    }

    static get none() {
        return new Stroke(null,null,null,false)
    }
}
class ShvgCanvas extends HTMLElement {
    constructor() {
        super();
        this.attachShadow({ mode: 'open' });

        this.svgRenderer = new ShvgRenderer(new Projection(new HyperbolicVector(0.001, 0.002),0,1,3.4,1));
        this.elements = [];
        this.shadowRoot.innerHTML = this.svgRenderer.generateBaseSVG();
    }

    point(vector, colour = "#000") {
        const point = this.svgRenderer.addPoint(vector, colour);
        this.elements.push(point);
    }

    geodesic(geo,stroke=new Stroke("#F00",5,true),samples=5) {
        if (stroke.enabled) {
            if (stroke.fixed) {
                this.elements.push(this.svgRenderer.addGeodesic(geo,stroke,samples));
            } else {
                let angle = geo.a.to(geo.b).angle
                let step = stroke.width/2
                let points = [geo.a.add(HyperbolicVector.fromPolar(step,angle+tau/4)),geo.a.add(HyperbolicVector.fromPolar(step,angle-tau/4)),geo.b.add(HyperbolicVector.fromPolar(step,angle-tau/4)),geo.b.add(HyperbolicVector.fromPolar(step,angle+tau/4))]
                this.elements.push(this.svgRenderer.addPolygon(points,samples,stroke.colour,Stroke.none));
            }
        }
    }
    
    polygon(points,samples=5, fill = "#0F0", stroke = new Stroke("#00F",5,true)) {
        this.elements.push(this.svgRenderer.addPolygon(points,samples,fill,stroke));
    }

    blob(points, fill = "#0F0", stroke = new Stroke("#00F",5,true)) {
        this.elements.push(this.svgRenderer.addBlob(points,fill,stroke));
    }

    circle(centre,radius=1,samples=10,fill="#888",stroke=new Stroke("#FFF",0.03,false)) {
        if (stroke.enabled)
        {
            if (stroke.fixed)
            {
                this.elements.push(this.svgRenderer.addCircle(centre,radius,samples,fill,stroke));
            } else {
                this.elements.push(this.svgRenderer.addCircle(centre,radius+stroke.width,samples,stroke.colour,Stroke.none));
                this.elements.push(this.svgRenderer.addCircle(centre,radius,samples,fill,stroke));
            }
        } else {
            this.elements.push(this.svgRenderer.addCircle(centre,radius,samples,fill,stroke));
        }
    }

    clear() {
        this.elements = [];
        this.shadowRoot.innerHTML = this.svgRenderer.generateBaseSVG();
    }

    draw() {
        const svgContent = this.svgRenderer.render(this.elements);
        this.shadowRoot.innerHTML = svgContent;
    }
}

customElements.define('shvg-canvas', ShvgCanvas);

class polygon {
    constructor(points) {
        this.points=points
    }
}

function triangleSideLength(AB,BC,AC) { //Based on the hyperbolic law of cosines
    return Math.acosh((Math.cos(BC)+Math.cos(AB)*Math.cos(AC))/(Math.sin(AB)*Math.sin(AC)))
}

function createTiling(p, q, maxDepth) {
    let points = []
    let edges = []

    const centreToVertex = triangleSideLength(tau/p,tau/(2*q),tau/(2*q))
    const centerToMidpoint = triangleSideLength(tau/(2*p),tau/(q*2),tau/4)

    let layer = [new polygon([])]

    let prev = HyperbolicVector.fromPolar(centreToVertex,(p-1)*(tau/p))
    for (let i=0;i<p;i++)
    {
        let p1 = HyperbolicVector.fromPolar(centreToVertex,i*(tau/p))
        points.push(p1)
        layer[0].points.push(p1)
        edges.push([p1,prev])
        prev = p1
    }

    function reflect(poly,direction,p) {
        let p1=poly.points[(direction+1)%p]
        let p2=poly.points[direction]

        let geo = new Geodesic(p1,p2)

        let others = []
        for (let i = 1; i<(p-1);i++) 
        {
            let index = (direction+p-i)
            let newPoint = poly.points[index%p].reflect(geo)
            others.push(newPoint)
            points.push(newPoint)
        }

        poly = new polygon([p1,p2].concat(others))

        poly.points.forEach((point,i) => {
            edges.push([point,poly.points[(i+1)%p]])
        });

        return poly
    }
    
    for (let depth = 0; depth<maxDepth; depth++)
    {
        let nextLayer = []
        
        layer.forEach(poly => {
            for (let i=0;i<p;i++)
            {
                nextLayer.push(reflect(poly,i,p))
            }
        })
        
        layer = nextLayer
    }

    return {points, edges}
}

// function createSubTiling(p, q, maxDepth,subRadius) {
//     let points = []
//     let edges = []

//     const centreToVertex = triangleSideLength(tau/p,tau/(2*q),tau/(2*q))
//     const centerToMidpoint = triangleSideLength(tau/(2*p),tau/(q*2),tau/4)

//     let layer = [new polygon([])]

//     let prev = HyperbolicVector.fromPolar(centreToVertex,(p-1)*(tau/p))
//     for (let i=0;i<p;i++)
//     {
//         let p1 = HyperbolicVector.fromPolar(centreToVertex,i*(tau/p))
//         points.push(p1.subject(subRadius))
//         layer[0].points.push(p1)
//         edges.push([p1.subject(subRadius),prev.subject(subRadius)])
//         prev = p1
//     }

//     function reflect(poly,direction,p) {
//         let p1=poly.points[(direction+1)%p]
//         let p2=poly.points[direction]

//         let geo = new Geodesic(p1,p2)

//         let others = []
//         for (let i = 1; i<(p-1);i++) 
//         {
//             let index = (direction+p-i)
//             let newPoint = poly.points[index%p].reflect(geo)
//             others.push(newPoint)
//             points.push(newPoint.subject(subRadius))
//         }

//         poly = new polygon([p1,p2].concat(others))

//         poly.points.forEach((point,i) => {
//             edges.push([point.subject(subRadius),poly.points[(i+1)%p].subject(subRadius)])
//         });

//         return poly
//     }
    
//     for (let depth = 0; depth<maxDepth; depth++)
//     {
//         let nextLayer = []
        
//         layer.forEach(poly => {
//             for (let i=0;i<p;i++)
//             {
//                 nextLayer.push(reflect(poly,i,p))
//             }
//         })
        
//         layer = nextLayer
//     }

//     reflect(layer[0],1,p)
//     reflect(layer[1],1,p)
//     reflect(layer[2],1,p)
//     reflect(layer[3],1,p)
//     reflect(layer[4],1,p)
//     reflect(layer[5],1,p)

//     return {points, edges}
// }
function createSubTiling(p, q, maxDepth,subRadius) {
    let points = []
    let edges = []

    const centreToVertex = triangleSideLength(tau/p,tau/(2*q),tau/(2*q))
    const centerToMidpoint = triangleSideLength(tau/(2*p),tau/(q*2),tau/4)

    let layer = [new polygon([])]

    let prev = HyperbolicVector.fromPolar(centreToVertex,(p-1)*(tau/p))
    for (let i=0;i<p;i++)
    {
        let p1 = HyperbolicVector.fromPolar(centreToVertex,i*(tau/p))
        points.push(p1.scale(subRadius))
        layer[0].points.push(p1)
        edges.push([p1.scale(subRadius),prev.scale(subRadius)])
        prev = p1
    }

    function reflect(poly,direction,p) {
        let p1=poly.points[(direction+1)%p]
        let p2=poly.points[direction]

        let geo = new Geodesic(p1,p2)

        let others = []
        for (let i = 1; i<(p-1);i++) 
        {
            let index = (direction+p-i)
            let newPoint = poly.points[index%p].reflect(geo)
            others.push(newPoint)
            points.push(newPoint.scale(subRadius))
        }

        poly = new polygon([p1,p2].concat(others))

        poly.points.forEach((point,i) => {
            edges.push([point.scale(subRadius),poly.points[(i+1)%p].scale(subRadius)])
        });

        return poly
    }
    
    for (let depth = 0; depth<maxDepth; depth++)
    {
        let nextLayer = []
        
        layer.forEach(poly => {
            for (let i=0;i<p;i++)
            {
                nextLayer.push(reflect(poly,i,p))
            }
        })
        
        layer = nextLayer
    }

    reflect(layer[0],1,p)
    reflect(layer[1],1,p)
    reflect(layer[2],1,p)
    reflect(layer[3],1,p)
    reflect(layer[4],1,p)
    reflect(layer[5],1,p)

    return {points, edges}
}
/*
complex fast fourier transform and inverse from
http://rosettacode.org/wiki/Fast_Fourier_transform#C.2B.2B
*/
function icfft(amplitudes) {
    var N = amplitudes.length;
    var iN = 1 / N;

    //conjugate if imaginary part is not 0
    for(var i = 0 ; i < N; ++i)
        if(amplitudes[i] instanceof Complex)
            amplitudes[i].im = -amplitudes[i].im;

    //apply fourier transform
    amplitudes = cfft(amplitudes)

    for(var i = 0 ; i < N; ++i)
    {
        //conjugate again
        amplitudes[i].im = -amplitudes[i].im;
        //scale
        amplitudes[i].re *= iN;
        amplitudes[i].im *= iN;
    }
    return amplitudes;
}

function cfft(amplitudes) {
    var N = amplitudes.length;
    if( N <= 1 )
        return amplitudes;

    var hN = N / 2;
    var even = [];
    var odd = [];
    even.length = hN;
    odd.length = hN;
    for(var i = 0; i < hN; ++i)
    {
        even[i] = amplitudes[i*2];
        odd[i] = amplitudes[i*2+1];
    }
    even = cfft(even);
    odd = cfft(odd);

    var a = -2*Math.PI;
    for(var k = 0; k < hN; ++k)
    {
        if(!(even[k] instanceof Complex))
            even[k] = new Complex(even[k], 0);
        if(!(odd[k] instanceof Complex))
            odd[k] = new Complex(odd[k], 0);
        var p = k/N;
        var t = new Complex(0, a * p);
        t.cexp(t).mul(odd[k], t);
        amplitudes[k] = even[k].add(t, odd[k]);
        amplitudes[k + hN] = even[k].sub(t, even[k]);
    }
    return amplitudes;
}


/*
basic complex number arithmetic from
http://rosettacode.org/wiki/Fast_Fourier_transform#Scala
*/
function Complex(re, im) {
    this.re = re;
    this.im = im || 0.0;
}
Complex.prototype.add = function(other, dst) {
    dst.re = this.re + other.re;
    dst.im = this.im + other.im;
    return dst;
}
Complex.prototype.sub = function(other, dst) {
    dst.re = this.re - other.re;
    dst.im = this.im - other.im;
    return dst;
}
Complex.prototype.mul = function(other, dst) {
    //cache re in case dst === this
    var r = this.re * other.re - this.im * other.im;
    dst.im = this.re * other.im + this.im * other.re;
    dst.re = r;
    return dst;
}
Complex.prototype.cexp = function(dst) {
    var er = Math.exp(this.re);
    dst.re = er * Math.cos(this.im);
    dst.im = er * Math.sin(this.im);
    return dst;
}
Complex.prototype.log = function() {
    /*
    although 'It's just a matter of separating out the real and imaginary parts of jw.' is not a helpful quote
    the actual formula I found here and the rest was just fiddling / testing and comparing with correct results.
    http://cboard.cprogramming.com/c-programming/89116-how-implement-complex-exponential-functions-c.html#post637921
    */
    if( !this.re )
        console.log(this.im.toString()+'j');
    else if( this.im < 0 )
        console.log(this.re.toString()+this.im.toString()+'j');
    else
        console.log(this.re.toString()+'+'+this.im.toString()+'j');
}


// Ajouté par valentin
// vvvvvvvvvvvvvvvvvvv

Complex.prototype.modSquared = function (){
    return this.re * this.re + this.im * this.im;
}

Complex.prototype.mod = function (){
    return math.sqrt(this.modSquared());
}




// Simplifications pour le projet trac


/**
 * Calcule la transformée de fourier de signal et la présente de façon à ce qu’elle soit facile à traiter (bonnes amplitudes à toutes les fréquences & suppression de la symétrie)
 * @param signal
 * @param fe féquence d’échantillonnage ( 1/(temps entre 2 échantillons) )
 * @returns un tableau a 2 dimensions : [ fréquences, amplitudes]
 */
function fft(signal, fe){

    // ajuster la taille du signal à une puissance de 2 (zéro padding)
    var puissanceDeDeux = Math.ceil(Math.log2(signal.length));
    var nbEchAjouter = Math.pow(2,puissanceDeDeux) - signal.length;

    var zeros = Array(nbEchAjouter).fill(0);
    var signalATraiter = [...signal].concat(zeros); // todo virer la copie de tableau ?

    var N = signalATraiter.length;
    var fft_cplx = cfft(signalATraiter);
    var range_f = Array(N / 2);
    var fft_abs = Array(N / 2);

    var pasFreq = fe/N;

    for(var i=0 ; i<(N / 2) ; i++){

        range_f[i] = pasFreq * i;
        if(i===0){
            fft_abs[i] = fft_cplx[i].mod()/N;
        }else{
            fft_abs[i] = fft_cplx[i].mod() * 2/N;
        }
    }

    return [range_f,fft_abs];

}

function detectionPic(signal, tailleTroncon, seuil) {
    var nbTroncons = Math.ceil(signal.length / tailleTroncon);
    var maxValeurTroncon = Array(nbTroncons).fill(0);
    var maxPositionTroncon = Array(nbTroncons).fill(0);
    var max=[];


    for (var t = 0; t < nbTroncons; t++) {
        for (var i = 0; (i < tailleTroncon) && ((t*tailleTroncon  + i) < signal.length); i++) {
            if (maxValeurTroncon[t] < signal[t*tailleTroncon  + i]) {
                maxValeurTroncon[t] = signal[t*tailleTroncon  + i];
                maxPositionTroncon[t] = t * tailleTroncon + i;
            }
        }

        if(t>=1){
            var t2;
            if(t>1){

                t2= maxValeurTroncon[t-2];
            }else {
                t2=0;
            }
            var t1= maxValeurTroncon[t-1];
            var t0= maxValeurTroncon[t];
            if(Math.abs(t2 - t1)>seuil && Math.abs(t1-t0)>seuil){
                max.push([t1,maxPositionTroncon[t-1]]);
            }
        }
    }

    return [maxPositionTroncon,maxValeurTroncon,max];
}

function traceTroncons(signalLength, max, tailleTroncon){

    var trace= Array(signalLength);
    var nbTroncons = Math.ceil(signalLength / tailleTroncon);

    for (var t = 0; t < nbTroncons; t++) {
        for (var i = 0; (i < tailleTroncon) && ((t*tailleTroncon  + i) < signalLength); i++) {
            trace[t*tailleTroncon  + i] = max[t];
        }
    }

    return trace;
}

class TraiteurAcceleration {

    constructor(nbPoints, te) {
        this.acceleration = [new BufferCirculaire(nbPoints), new BufferCirculaire(nbPoints), new BufferCirculaire(nbPoints), new BufferCirculaire(nbPoints)]; // x, y, z, norme
        this.vitesse = [new BufferCirculaire(nbPoints), new BufferCirculaire(nbPoints), new BufferCirculaire(nbPoints), new BufferCirculaire(nbPoints)];      //idem
        this.position = [new BufferCirculaire(nbPoints), new BufferCirculaire(nbPoints), new BufferCirculaire(nbPoints), new BufferCirculaire(nbPoints)];     //idem

        this.buffersLength = nbPoints;

        this.te = te;
    }

    ajouterEchantillon(x, y, z){
        this.acceleration[0].push(x);
        this.acceleration[1].push(y);
        this.acceleration[2].push(z);

        var norme = Math.sqrt( Math.pow(x,2) + Math.pow(y,2) + Math.pow(z,2) );
        this.acceleration[3].push(norme);

        this.#miseAJourVitesse();
        this.#miseAJourPosition();

    }


    #miseAJour(buffer_data_array, buffer_derivee_array){
        for( var i = 0 ; i < 4 ; i++ ){
            var donnee = this.te * buffer_data_array[i].get(0) + buffer_derivee_array[i].get(0);//Attention vitesse pas encore mise à jour donc vitesse[i].get(0) => echantillion précédant
            buffer_derivee_array[i].push(donnee);
        }
    }
    #miseAJourVitesse(){
        this.#miseAJour(this.acceleration, this.vitesse);
    }

    #miseAJourPosition() {
        this.#miseAJour(this.vitesse, this.position);
    }


    static #calculMoyenne(buffer, nbPoints){
        if( typeof(nbPoints) == 'undefined' ) {
            nbPoints = buffer.getSize();
        }
        var moyenne = 0;


        for (var i = 0; i < nbPoints; i++) {
            moyenne += buffer.get(i);
        }

        moyenne /= nbPoints;

        return moyenne;
    }

    #accelerationMoyenne(composante, nbPoints){

        return TraiteurAcceleration.#calculMoyenne(this.acceleration[composante],nbPoints);
    }

    #vitesseMoyenne(composante, nbPoints){
        return TraiteurAcceleration.#calculMoyenne(this.vitesse[composante],nbPoints);
    }


    #positionMoyenne(composante, nbPoints){
        return TraiteurAcceleration.#calculMoyenne(this.position[composante],nbPoints);
    }



    /**
     *
     * @param type  0=acceleration, 1= vitesse, 2=position
     * @returns vecteur de tempp
     */
    #getTemps(type){
        var tailleTemps;
        switch (type) {
            case 0:
                tailleTemps = this.acceleration[0].getSize();
                break;
            case 1:
                tailleTemps = this.vitesse[0].getSize();
                break;
            case 2:
                tailleTemps = this.position[0].getSize();
                break;
        }

        var vecteurTemps = Array(tailleTemps);
        for( var i = 0 ; i < vecteurTemps.length ; i++){
            vecteurTemps[i] = this.te *(i - tailleTemps + 1 );
        }

        return vecteurTemps;

    }


    #fftAcceleration(composante, nbPoints) {
        return fft(this.acceleration[composante].enTableau(nbPoints), 1 / this.te);
    }



    #fftVitesse(composante, nbPoints) {
        return fft(this.vitesse[composante].enTableau(nbPoints), 1 / this.te);
    }

    #fftPosition(composante, nbPoints) {
        return fft(this.position[composante].enTableau(nbPoints), 1 / this.te);
    }

    static #str2composante(composante_str){
        switch (composante_str.toLowerCase()){
            case "x":
                return 0;
            case "y":
                return 1;
            case "z":
                return 2;
            case "norme":
            case "n":
                return 3;
            default:
                throw "composante inconnue peut etre «x»,«y»,«z» ou «norme»";
        }
    }

    getAcceleration(composanteStr){
        return this.acceleration[ TraiteurAcceleration.#str2composante(composanteStr) ];
    }

    getVitesse(composanteStr){
        return this.vitesse[ TraiteurAcceleration.#str2composante(composanteStr) ];
    }

    getPosition(composanteStr){
        return this.position[ TraiteurAcceleration.#str2composante(composanteStr) ];
    }

    getTempsAcceleration() {
        return this.#getTemps(0);
    }
    getTempsVitesse() {
        return this.#getTemps(1);
    }
    getTempsPosition() {
        return this.#getTemps(2);
    }

    getAccelerationMoyenne(composanteStr, nbPoints){
        return this.#accelerationMoyenne( TraiteurAcceleration.#str2composante(composanteStr), nbPoints );
    }

    getVitesseMoyenne(composanteStr, nbPoints){
        return this.#vitesseMoyenne( TraiteurAcceleration.#str2composante(composanteStr), nbPoints );
    }

    getPositionMoyenne(composanteStr, nbPoints){
        return this.#positionMoyenne( TraiteurAcceleration.#str2composante(composanteStr), nbPoints );
    }

    getfftAcceleration(composanteStr,nbPoints){
        return this.#fftAcceleration( TraiteurAcceleration.#str2composante(composanteStr), nbPoints );
    }

    getfftVitesse(composanteStr){
        return this.#fftVitesse( TraiteurAcceleration.#str2composante(composanteStr) );
    }

    getfftPosition(composanteStr){
        return this.#fftPosition( TraiteurAcceleration.#str2composante(composanteStr) );
    }
}


class BufferCirculaire {

    constructor(nbElts) {
        this.buffer = Array(nbElts).fill(0);
        this.size = 0;
        this.index = nbElts - 1;
    }

    getSize() {
        return this.size;
    }

    get(n){
        return this.buffer[Math.abs(n - this.index - this.buffer.length) % this.buffer.length];
    }

    push(value){
        this.index = Math.abs((this.index + 1) % this.buffer.length);
        this.buffer[this.index] = value;
        if(this.size<this.buffer.length){
            this.size++;
        }
    }


    enTableau(nbPoints){
        if( typeof(nbPoints) == 'undefined' ){
            nbPoints = this.getSize();
        }

        var tableau = Array(nbPoints);
        for( var i = 0 ; i < tableau.length ; i++ ){
            tableau[tableau.length-1-i] = this.get(i);
        }

        return(tableau);
    }



}
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

// todo renommer en spectre? ( fft trompeur maintenant ? …)
/**
 * Calcule la transformée de fourier de signal et la présente de façon à ce qu’elle soit
 * facile à traiter (bonnes amplitudes à toutes les fréquences & suppression de la symétrie)
 *
 * @param signal enregistrement à traiter
 *
 * @param fe féquence d’échantillonnage ( 1/(temps entre 2 échantillons) )
 *
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


/**
 * Découpe l’enregistrement en tronçons et trouve les maximums de chaque tronçon puis pour chaque tronçon,
 * on compare le maximum de chaque tronçon et si sa différence par rapport à ses voisins est importante,
 * on ajoute son maximum et la location du maximum dans la liste des maximums
 *
 * @param signal enregistrement à traiter
 *
 * @param tailleTroncon taille d’un tronçon en nombre d’échantillons
 *
 * @param seuil différence minimale entre un tronçon et ses voisins pour qu’il soit considèré comme un maximum
 *
 * @returns [
 *              tableau la position des maximums pour chaque tronçon,
 *              tableau contenant les maximums de chaque tronçon,
 *              [
 *                  [pos1,max1],
 *                  [pos2,max2],
 *                  …
 *                  [posN,maxN],
 *              ]
 *          ]
 */
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


/**
 * génère un vecteur de points pour tracer les maximums (utilisé pour les tests)
 *
 * @param signalLength taille de l’enregistrement passé à la fonction detectionPic
 *
 * @param max tableau retourné par detectionPic
 *
 * @param tailleTroncon taille d’un tronçon passé à la fonction detectionPic
 *
 * @returns vecteur de points
 */
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

/**
 * Classe utilisée pour le traitement des données de l’acceleromètre
 */
class TraiteurAcceleration {


    /**
     * @param nbPoints nombre de points à enregistrer
     *
     * @param te periode d’échantillonnage
     */
    constructor(nbPoints, te) {
        this.acceleration = [new BufferCirculaire(nbPoints), new BufferCirculaire(nbPoints), new BufferCirculaire(nbPoints), new BufferCirculaire(nbPoints)]; // x, y, z, norme
        this.vitesse = [new BufferCirculaire(nbPoints), new BufferCirculaire(nbPoints), new BufferCirculaire(nbPoints), new BufferCirculaire(nbPoints)];      //idem
        this.position = [new BufferCirculaire(nbPoints), new BufferCirculaire(nbPoints), new BufferCirculaire(nbPoints), new BufferCirculaire(nbPoints)];     //idem

        this.buffersLength = nbPoints;

        this.te = te;
    }


    /**
     * Ajoute un échantillon récupere par l’acceléromètre aux données enregistrées et met à jour la vitesse et la position
     * @param x
     * @param y
     * @param z
     */
    ajouterEchantillon(x, y, z){
        this.acceleration[0].push(x);
        this.acceleration[1].push(y);
        this.acceleration[2].push(z);

        var norme = Math.sqrt( Math.pow(x,2) + Math.pow(y,2) + Math.pow(z,2) );
        this.acceleration[3].push(norme);

        this.#miseAJourVitesse();
        this.#miseAJourPosition();

    }


    /**
     * Met à jour les buffers de buffer_intg_array en integrant le contenu de buffer_src_array
     *
     * @param buffer_src_array
     *
     * @param buffer_intg_array
     */
    #miseAJour(buffer_src_array, buffer_intg_array){
        for( var i = 0 ; i < 4 ; i++ ){
            var donnee = this.te * buffer_src_array[i].get(0) + buffer_intg_array[i].get(0);//Attention vitesse pas encore mise à jour donc vitesse[i].get(0) => echantillion précédant
            buffer_intg_array[i].push(donnee);
        }
    }

    /**
     * Met à jour la vitesse à partir du dernier échantillion de l’accélération ajouté
     */
    #miseAJourVitesse(){
        this.#miseAJour(this.acceleration, this.vitesse);
    }

    /**
     * Met à jour la position à partir du dernier échantillion de vitesse ajouté
     *
     * <em> la fonction miseAJourVitesse doit etre appelée avant</em>
     */
    #miseAJourPosition() {
        this.#miseAJour(this.vitesse, this.position);
    }


    /**
     * Calcule la valeur moyenne des nbPoints derniers points du buffer buffer
     *
     * @param buffer signal à traiter
     *
     * @param nbPoints nombre d’échantillons à traiter
     *
     * @returns la valeur moyenne du signal sur nbPoints points
     */
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

    /**
     * Calcule l’accélération moyenne sur nbPoints points sur la composante composante
     *
     * @param composante 0=x, 1=y, 2=z, 3=norme
     *
     * @param nbPoints nombre de points
     *
     * @returns valeur moyenne de l’acceleration sur l’axe choisi
     */
    #accelerationMoyenne(composante, nbPoints){
        return TraiteurAcceleration.#calculMoyenne(this.acceleration[composante],nbPoints);
    }

    /**
     * Calcule la vitesse moyenne sur nbPoints points sur la composante composante
     *
     * @param composante 0=x, 1=y, 2=z, 3=norme
     *
     * @param nbPoints nombre de points
     *
     * @returns valeur moyenne de la vitesse sur l’axe choisi
     */
    #vitesseMoyenne(composante, nbPoints){
        return TraiteurAcceleration.#calculMoyenne(this.vitesse[composante],nbPoints);
    }

    /**
     * Calcule la position moyenne sur nbPoints points sur la composante composante
     *
     * @param composante 0=x, 1=y, 2=z, 3=norme
     *
     * @param nbPoints nombre de points
     *
     * @returns valeur moyenne de la position sur l’axe choisi
     */
    #positionMoyenne(composante, nbPoints){
        return TraiteurAcceleration.#calculMoyenne(this.position[composante],nbPoints);
    }



    /**
     * Un vecteur temps correspondant au type de données choisies
     *
     * @param type  0=acceleration, 1= vitesse, 2=position
     *
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


    /**
     * Calcule le spectre de l’accélération sur la compsante composante
     *
     * @param composante 0=x, 1=y, 2=z, 3=norme
     *
     * @param nbPoints nombre de points sur lesquels calculer la fft ( de préférence une puissance de 2)
     *
     * @returns un tableau a 2 dimensions : [ fréquences, amplitudes]
     */
    #fftAcceleration(composante, nbPoints) {
        return fft(this.acceleration[composante].enTableau(nbPoints), 1 / this.te);
    }


    /**
     * Calcule le spectre de la vitesse sur la compsante composante
     *
     * @param composante 0=x, 1=y, 2=z, 3=norme
     *
     * @param nbPoints nombre de points sur lesquels calculer la fft ( de préférence une puissance de 2)
     *
     * @returns un tableau a 2 dimensions : [ fréquences, amplitudes]
     */
    #fftVitesse(composante, nbPoints) {
        return fft(this.vitesse[composante].enTableau(nbPoints), 1 / this.te);
    }

    /**
     * Calcule le spectre de la position sur la compsante composante
     *
     * @param composante 0=x, 1=y, 2=z, 3=norme
     *
     * @param nbPoints nombre de points sur lesquels calculer la fft ( de préférence une puissance de 2)
     *
     * @returns un tableau a 2 dimensions : [ fréquences, amplitudes]
     */
    #fftPosition(composante, nbPoints) {
        return fft(this.position[composante].enTableau(nbPoints), 1 / this.te);
    }

    /**
     * Convertit la chaine de carractère composante_str en index pour les tableaux this.acceleration, this.vitesse, this.position
     *
     * @param composante_str
     *
     * @returns index
     */
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

    /**
     * @param composanteStr
     *
     * @returns Le buffer acceleration correspondant à la composante composante
     */
    getAcceleration(composanteStr){
        return this.acceleration[ TraiteurAcceleration.#str2composante(composanteStr) ];
    }

    /**
     * @param composanteStr
     *
     * @returns Le buffer vitesse correspondant à la composante composante
     */
    getVitesse(composanteStr){
        return this.vitesse[ TraiteurAcceleration.#str2composante(composanteStr) ];
    }

    /**
     * @param composanteStr
     *
     * @returns Le buffer position correspondant à la composante composante
     */
    getPosition(composanteStr){
        return this.position[ TraiteurAcceleration.#str2composante(composanteStr) ];
    }

    /**
     * @returns un vecteur temps correspondant à l’acceleration
     */
    getTempsAcceleration() {
        return this.#getTemps(0);
    }

    /**
     * @returns un vecteur temps correspondant à la vitesse
     */
    getTempsVitesse() {
        return this.#getTemps(1);
    }

    /**
     * @returns un vecteur temps correspondant à la position
     */
    getTempsPosition() {
        return this.#getTemps(2);
    }

    /**
     * Retourne l’acceleration moyenne de l’acceleration
     *
     * @param composanteStr composante
     *
     * @param nbPoints nombre de points sur lequels calculer l’acceleration moyenne
     *
     * @returns l’acceleration moyenne de la composant composanteStr sur nbPoints points
     */
    getAccelerationMoyenne(composanteStr, nbPoints){
        return this.#accelerationMoyenne( TraiteurAcceleration.#str2composante(composanteStr), nbPoints );
    }

    /**
     * Retourne la vitesse moyenne
     *
     * @param composanteStr composante
     *
     * @param nbPoints nombre de points sur lequels calculer la vitesse moyenne
     *
     * @returns la vitesse moyenne de la composant composanteStr sur nbPoints points
     */
    getVitesseMoyenne(composanteStr, nbPoints){
        return this.#vitesseMoyenne( TraiteurAcceleration.#str2composante(composanteStr), nbPoints );
    }

    /**
     * Retourne la position moyenne
     *
     * @param composanteStr composante
     *
     * @param nbPoints nombre de points sur lequels calculer la vitesse moyenne
     *
     * @returns la vitesse moyenne de la composant composanteStr sur nbPoints points
     */
    getPositionMoyenne(composanteStr, nbPoints){
        return this.#positionMoyenne( TraiteurAcceleration.#str2composante(composanteStr), nbPoints );
    }

    /**
     * Retourne le spectre de l’accélération
     *
     * @param composanteStr composante
     *
     * @param nbPoints nombre de points sur lequels calculer la fft
     *
     * @returns un tableau a 2 dimensions : [ fréquences, amplitudes]
     */
    getfftAcceleration(composanteStr,nbPoints){
        return this.#fftAcceleration( TraiteurAcceleration.#str2composante(composanteStr), nbPoints );
    }

    /**
     * Retourne le spectre de la vitesse
     *
     * @param composanteStr composante
     *
     * @param nbPoints nombre de points sur lequels calculer la fft
     *
     * @returns un tableau a 2 dimensions : [ fréquences, amplitudes]
     */
    getfftVitesse(composanteStr, nbPoints){
        return this.#fftVitesse( TraiteurAcceleration.#str2composante(composanteStr), nbPoints);
    }

    /**
     * Retourne le spectre de la position
     *
     * @param composanteStr composante
     *
     * @param nbPoints nombre de points sur lequels calculer la fft
     *
     * @returns un tableau a 2 dimensions : [ fréquences, amplitudes]
     */
    getfftPosition(composanteStr, nbPoints){
        return this.#fftPosition( TraiteurAcceleration.#str2composante(composanteStr), nbPoints );
    }
}


class BufferCirculaire {

    constructor(nbElts) {
        this.buffer = Array(nbElts).fill(0);
        this.size = 0;
        this.index = nbElts - 1;
    }

    /**
     * Retourne le nombre d’éléments mis dans le buffer circulaire
     * @returns {number}
     */
    getSize() {
        return this.size;
    }

    /**
     *
     * @param n
     * @returns le n ième dernier élément ajouté dans le tableau
     */
    get(n){
        return this.buffer[Math.abs(n - this.index - this.buffer.length) % this.buffer.length];
    }

    /**
     * Ajoute value au buffer circulaire si la taille le buffer est plein, efface la valeur la plus ancienne
     *
     * @param value
     */
    push(value){
        this.index = Math.abs((this.index + 1) % this.buffer.length);
        this.buffer[this.index] = value;
        if(this.size<this.buffer.length){
            this.size++;
        }
    }


    /**
     * Retourne un tableau contenant les valeurs ajoutées au buffer
     *
     * @param nbPoints taille du tableau
     *
     * @returns un tableau des nbPoints dernières valeurs ajoutées
     */
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
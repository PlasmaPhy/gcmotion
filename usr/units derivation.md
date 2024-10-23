---
title: Λαγκραντζιανός και Χαμιλτονιανός Φορμαλισμός Guiding Center
author: George Tsiamasiotis
date: Last compiled on \\today

geometry: left=20mm,right=20mm,top=20mm,bottom=25mm
toc: true
fontsize: 10pt
lang: gr
mainfont: Kerkis
latex-engine: xelatex
---

## Προαπαιτούμενα

### Περί συντεταγμένων

Η επιλογή των συντεταγμένων γενικά εξαρτάται άμεσα από τη μορφή του εκάστοτε μαγνητικού πεδίου. Τα συστήματα συντεταγμένων αυτά, παρότι δυσνόητα, διευκολύνουν σημαντικά την ανάλυση. Παρόλα αυτά, στο τέλος πρέπει να μετατραπούν σε συντεταγμένες εργαστηρίου.

Χρησιμοποιούμε το γνωστό σύστημα $(\psi,\theta,\zeta)$. Η συντεταγμένη $\psi$ ονομάζεται *τοροειδής ροή* (*toroidal flux*) ορίζει τις εμφωλευμένες (nested) καμπύλες, οι οποίες ορίζονται μέσω της $\psi = \text{σταθ}$. Το $\psi$ **πρέπει** να είναι $0$ στον μαγνητικό άξονα, και να αυξάνεται προς τα έξω. Άρα η $\psi$ είναι μη αρνητική.

Η μεταβλητή $\psi_p$ ονομάζεται *πολοειδής ροή* (*poloidal flux*) και σχετίζεται με την $\psi$ μέσω του safety factor q:

$$
q(\psi) = \dfrac{d\psi}{d\psi_p}
$$

Οι δύο μεταβλητές $\psi,\psi_p$ μπορούν να χρησιμοποιηθούν εναλλάξ, ανάλογα με το πρόβλημα.

Η γενική μορφή του μαγνητικού πεδίου είναι:

$$
\vec B = \nabla\psi \times \nabla\theta - \nabla\psi_p \times \nabla\zeta
$$

Εδώ φαίνεται ο Χαμιλτονιανός χαρακτήρας του μαγνητικού πεδίου, καθώς οι μαγνητικές γραμμές διέπονται από τις κανονικές εξισώσεις Hamilton:

$$
\dfrac{d\psi}{d\zeta} = -\partial_\theta\psi_p,\qquad\text{και}\qquad 
\dfrac{d\psi_p}{d\zeta} = \partial_\psi\psi_p
$$

Δηλαδή η $\psi_p(\psi,\theta,\zeta)$ είναι η Χαμιλτονιανή, $\psi$ η κανονική ορμή, $\theta$ η θέση και $\zeta$ ο χρόνος.

#### White-Boozer Coordinates

Εμείς χρησιμοποιούμε τις λεγόμενες *White-Boozer Coordinates*, στις οποίες οι μαγνητικές επιφάνειες εκφράζονται ως:

- Ανταλλοίωτη (contravariant) μορφή:  $\qquad \vec B = \nabla\psi \times \nabla\theta - \nabla\psi_p \times \nabla\zeta$

- Συναλλοίωτη (covariant) μορφή: $\qquad\quad\ \  \vec B = g(\psi)\nabla\zeta + I(\psi)\nabla\theta + \delta(\psi,\theta)\nabla\psi$

όπου:

- $g(\psi)$ το (μέσο) πολοειδές ρεύμα στο **εξωτερικό** της καμπύλης $\psi$,

- $I(\psi)$ το (μέσο) τοροειδές ρεύμα στο **εσωτερικό** της καμπύλης $\psi$,

- $\delta$ ένας όρος που συνδέεται με την ορθογωνιότητα του συστήματος συντεταγμένων, και εμείς θεωρούμε 0.

***

Για το μαγνητικό πεδίο ισχύει ο νόμος του Ampere, ο οποίος στο SI γράφεται:

$$
\vec\nabla\times\vec B = \mu_0\vec J
$$

Επομένως, αν στο SI το $g$ έχει πράγματι διαστάσεις ρεύματος, η ποσότητα $g\nabla\zeta$ πρέπει να έχει διαστάσεις $[A/m]$ και άρα το γινόμενο $\mu_0g\nabla\zeta$ έχει μονάδες $[T]$.

Άρα η συναλλοίωτη μορφή του μαγνητικού πεδίου στο SI γραφεται:

- Συναλλοίωτη μορφή (SI): $\qquad\quad\ \  \vec B = \mu_0\bigg[g(\psi)\nabla\zeta + I(\psi)\nabla\theta + \delta(\psi,\theta)\nabla\psi\bigg]$

### Περί μονάδων

Σε όλο το βιβλίο χρησιμοποιούνται μονάδες κανονικοποιημένες ως προς το χαρακτηριστικό μήκος του εκάστοτε προβλήματος, και ως προς την τιμή του μαγνητικού πεδίου στον άξονα, ενώ παντού θεωρούμε $c=1$ (White σελ 5). Έτσι, και ο χρόνος κανονικοποιείται ως προς τον λεγόμετο *Alfven time*, ο οποίος ορίζεται τυπικά ως:

$$
\tau_A = \dfrac{\rho}{\upsilon_A}
$$

όπου $\rho$ η μικρή ακτίνα του τόρου και $\upsilon_A$ η *ταχύτητα Alfven*, δηλαδή η ταχύτητα διάδοσης των ηλεκτρομανγητικών κυμάτων στον όγκο του πλάσματος (Wikipedia).

## Εξαγωγή Χαμιλτονιανής (με διαστάσεις)

Ξεκινάμε από την γνωστές εκφράσεις για τη Λαγκραντζιανή και τη Χαμιλτονιανή ενός φορτισμένου σωματιδίου σε ηλεκτρομαγνητικό πεδίο, οι οποίες προκειμένου να δίνουν σωστά τον νόμο του Νεύτωνα μέσα από τις εξισώσεις E-L και τις κανονικές εξισώσεις Hamilton αντίστοιχα, πρέπει να έχουν τη μορφή:

$$
\large\boxed{
    \mathcal{L} = \bigg[e\vec A(\vec x,t) +\dfrac{1}{2}m \dot{\vec x}\bigg]\dot{\vec x} -e\Phi(\vec x,t),\qquad\text{και}\qquad
    \mathcal{H} = \dfrac{1}{2}m\dot{\vec x^2} + e\Phi(\vec x,t)
}
$$

Η μόνη υπόθεση που κάνουμε για το Β είναι ότι μεταβάλλεται **αργά** στο χώρο και στον χρόνο, σε σχέση με τις κλίμακες που ορίζονται από την γυροακτίνα και την κυκλοτρονική συχνότητα του σωματιδίου.

Συνδυάζοντας τις δύο παραπάνω, παίρνουμε την εξής έκφραση για τη Λαγκραντζιανή:

$$
\large\boxed{
    \mathcal{L} = \bigg[e\vec A(\vec x,t) +m \dot{\vec x}\bigg]\dot{\vec x} -\mathcal{H}(\dot{\vec x}, \vec x)
}
$$

Εξ' ορισμού, οι κανονικές ορμές θα είναι:

$$
p_i = \dfrac{\partial \cal L}{\partial \dot x_i} \Rightarrow
\boxed{p_i = m\dot x_i + qA_i}
$$

 - Οι κανονικές ορμές ***δεν*** είναι μετρήσιμες ποσότητες και δεν έχουν ***κανένα απολύτως φυσικό νόημα***. Είναι όμως χρήσιμες στους υπολογισμούς.

 - Η κινητική ορμή είναι η $P_i = p_i -qA_i=m\dot x_i$ και είναι μετρήσιμη ποσότητα.

 - Η ποσότητα $e\vec A$ λέγεται και *ορμή πεδίου*.

 - Οι εξισώσεις κίνησης θα δώσουν φυσικά $\vec \upsilon = \dot{\vec x}$. Η συνθήκη αυτή "επιλέγεται" από την μεταβολική αρχή. Παρόλο που θα ισχύει για την πραγματική κίνηση, δεν ισχύει για όλες τις δυνατές διαδρομές στον φασικό χώρο.

Θέτοντας λοιπόν την ταχύτητα του σωματιδίου ως $\vec \upsilon$, και χωρίς να ξεχνάμε ότι $\vec \upsilon = \dot{\vec x}$, μπορούμε να σπάσουμε την ταχύτητα σε δύο συνιστώσες, μία παράλληλη και μία κάθετη στη κατεύθυνση του πεδίου στη θέση εκείνη:

$$
\vec\upsilon = \upsilon_{||}\hat{b}+ w\hat c
$$

Ορίζουμε το διάνυσμα θέσης του guiding center $\vec X$ ως:

$$
\vec x = \vec X + \dfrac{mw}{|e|B}\hat a
$$

---

**ΣΗΜΕΙΩΣΗ**

Η ποσότητα $|e|$ έχει μονάδες $[C]$. Η απόλυτη τιμή σημαίνει ότι φεύγει μόνο το πρόσημο.

---

όπου $\hat a = \cos\xi \hat e_1 + \sin\xi \hat e_2$ το μοναδιαίο διάνυσμα θέσης του σωματιδίου ως προς το γυρόκεντρό του, και $mw/|e|B$ η ακτίνα Larmor του σωματιδίου. Στο σχήμα φαίνεται η παραπάνω γεωμετρία, υποθέτωντας ότι $m/|e| = 1$.

![Συντεταγμένες θέσης Guiding Center](./assets/images/GC_coordinates.png "Συντεταγμένες θέσης Guiding Center")

Τελικά η Λαγκραντζιανή γίνεται:

$$
\mathcal{L} = \bigg[e\vec A(\vec x,t) +m (\upsilon_{||}\hat{b}+ w\hat c)\bigg]
\bigg[\dot{\vec X} + \dfrac{d}{dt}\bigg( \dfrac{mw}{|e|B}\hat a \bigg) \bigg] -\mathcal{H}(\vec u, \vec X)
$$

Αναπτύσσοντας κατά Taylor την κάθε συνιστώσα του $\vec A$ ως προς $\vec x$ γύρω από το $\vec X$, και έχοντας κατά νου ότι 
$\vec\nabla\vec A = \dfrac{\partial A_1}{\partial x_1}+\dfrac{\partial A_2}{\partial x_2}+\dfrac{\partial A_3}{\partial x_3}$, έχουμε:

$$
\vec A(\vec x,t) \simeq \vec A(\vec X,t) +\overbrace{\bigg(\dfrac{mw}{|e|B}\hat a\bigg)}^{\vec x - \vec X}\vec\nabla\vec A(\vec X,t)
$$

Ξαναγράφουμε την Λαγκραντζιανή με $\vec A = \vec A(\vec X,t)$:

$$
\begin{aligned}
    \mathcal{L} &= \bigg[ e\bigg(\vec A + \dfrac{mw}{|e|B}\hat a\vec\nabla\vec A\bigg)
    +m (\upsilon_{||}\hat{b}+ w\hat c)\bigg]\bigg[\dot{\vec X} + \dfrac{d}{dt}\bigg( \dfrac{mw}{|e|B}\hat a \bigg) \bigg] -\mathcal{H}\\
    &=\bigg[ e\vec A + m(\upsilon_{||}\hat{b}+ w\hat c) \bigg]\dot{\vec X} 
    +e\dfrac{mw}{|e|B}\hat a\vec\nabla\vec A\bigg[ \dot{\vec X} + \dfrac{d}{dt}\bigg( \dfrac{mw}{|e|B}\hat a \bigg) \bigg] + 
    e\vec A\dfrac{d}{dt}\bigg( \dfrac{mw}{|e|B}\hat a \bigg) \\
    &\qquad\qquad\qquad\qquad\qquad\quad
    +m\upsilon_{||}\hat b \dfrac{d}{dt}\bigg( \dfrac{mw}{|e|B}\hat a \bigg) +
    mw\hat c \dfrac{d}{dt}\bigg( \dfrac{mw}{|e|B}\hat a \bigg) - \mathcal{H} \\
    &=\bigg[ e\vec A +m(\upsilon_{||}\hat{b}+ w\hat c)\bigg]\dot{\vec X} +
    \overbrace{e\vec A\dfrac{d}{dt}\bigg( \dfrac{mw}{|e|B}\hat a \bigg) +
    e\bigg( \dfrac{mw}{|e|B}\hat a \vec\nabla\bigg)\vec A \dot{\vec X}}^{\text{θα φύγουν με την }dS/dt} +
    e\bigg( \dfrac{mw}{|e|B}\hat a \vec\nabla\bigg)\vec A \dfrac{d}{dt}\bigg(\dfrac{mw}{|e|B}\hat a\bigg) \\
    &\qquad\qquad\qquad\qquad\quad\quad\quad
    +\overbrace{m\upsilon_{||}\hat b \dfrac{d}{dt}\bigg(\dfrac{mw}{|e|B}\hat a\bigg)}^{\text{0, αφού }d\hat a/dt \simeq\dot\xi\hat c = \mathcal{O}(\epsilon^2) \text{ και } \hat b \cdot\hat c = 0} 
        + m w\hat c\dfrac{d}{dt}\bigg(\dfrac{mw}{|e|B}\hat a\bigg) - \mathcal{H}\\
\end{aligned}
$$

---

**ΣΗΜΕΙΩΣΗ**

Εδώ υποθέτω ότι έτσι όπως έχει οριστεί η ακτίνα Larmor $r_L = \dfrac{mw}{|e|B}$ και η φορά του $\hat a$, η απόλυτη τιμή του $e$ φέυγει όταν η $r_l$ αποκοπεί από το $\hat a$. Συγκεκριμένα, έτσι όπως είναι ορισμένη η ακτίνα Larmor, δεν λαμβάνει υπ'οψιν το πρόσημο του σωματιδίου, και άρα η δουλειά πέφτει στο $\hat a$, δηλαδή πρέπει να έχει αντίθετη κατεύθυνση για ετερόσημα σωματίδια. Είναι εύλογο, αλλά μπορεί να είναι λάθος.

---

---

Ορίζουμε τη συνάρτηση $S = -e\dfrac{mw}{|e|B}\hat a\cdot\vec A$ και υπολογίζουμε την $\dfrac{dS}{dt}$:

$$
\begin{aligned}
    \dfrac{dS}{dt} &= -e\dfrac{dS}{dt} \bigg[ \dfrac{mw(\vec X,t)}{|e|B(\vec X,t)}\hat a(\vec X,t)\vec A(\vec X,t)\bigg]
    = -e\vec A \dfrac{d}{dt}\bigg(\dfrac{mw}{|e|B}\hat a\bigg) - e\bigg(\dfrac{mw}{|e|B}\hat a\bigg)\dfrac{d\vec A}{dt} = \\
    &=-e\vec A \dfrac{d}{dt}\bigg(\dfrac{mw}{|e|B}\hat a\bigg) - e\bigg(\dfrac{mw}{|e|B}\hat a\bigg)\bigg(\vec\nabla \vec A \dot{\vec X} + \dfrac{\partial \vec A}{dt}\bigg) =\\
    &=\overbrace{-e\vec A \dfrac{d}{dt}\bigg(\dfrac{mw}{|e|B}\hat a\bigg)- e\bigg(\dfrac{mw}{|e|B}\hat a\bigg)\vec\nabla \vec A \dot{\vec X}}^{\text{αυτά θα φύγουν με την }\mathcal{L}} - e\dfrac{mw}{|e|B}\hat a\cdot\dfrac{\partial \vec A}{dt} \\
\end{aligned}
$$

Υπολογίζουμε την νέα Λαγκραντζιανή: $\mathcal{L} = \mathcal{L}_{\text{παλιά}} + \dfrac{dS}{dt}$

$$
\begin{aligned}
    \mathcal{L} &= \bigg[ e\vec A +m(\upsilon_{||}\hat{b}+ w\hat c)\bigg]\dot{\vec X} 
    +e\bigg( \dfrac{mw}{|e|B}\hat a \vec\nabla\bigg)\vec A \dfrac{d}{dt}\bigg(\dfrac{mw}{|e|B}\hat a\bigg) 
    + m w\hat c\dfrac{d}{dt}\bigg(\dfrac{mw}{|e|B}\hat a\bigg)
    - e\dfrac{mw}{|e|B}\hat a\cdot\dfrac{\partial \vec A}{dt}\\
\end{aligned} 
$$

Yπολογίζουμε τους όρους ξεχωριστά, και χρησιμοποιούμε το γεγονός ότι $\dfrac{e}{|e|^2} = \dfrac{1}{e}$

$$
\begin{aligned}
    \bullet\quad
    mw\hat c\dfrac{d}{dt}\bigg(\dfrac{mw}{|e|B}\hat a\bigg) &= 
    \dfrac{m}{|e|}mw\hat c \bigg[ \dfrac{d}{dt}\bigg(\dfrac{w}{B}\bigg)\hat a + \dfrac{w}{B}\dfrac{d\hat a}{dt} \bigg] 
    \stackrel{d\hat a/dt = \dot\xi\hat c}{=} 
    =mw\dfrac{mw}{|e|B}\hat c\cdot\dot\xi\hat c  \\
    &\stackrel{\text{ΣΗΜΕΙΩΣΗ}}{=}mw\dfrac{mw}{qB}\dot\xi
    =\dfrac{m^2}{e}\dfrac{w^2}{B}\dot\xi \\
\end{aligned}
$$

> Η ποσότητα αυτή έχει μονάδες $[kgm^2/s^2]$, δηλαδή ενέργειας, όπως θα έπρεπε.

$$
\begin{aligned}
    \bullet\quad
    e\bigg( \dfrac{mw}{|e|B}\hat a \vec\nabla\bigg)\vec A \dfrac{d}{dt}\bigg(\dfrac{mw}{|e|B}\hat a\bigg) 
    &= e\dfrac{mw}{|e|B}(\hat a\vec\nabla)\vec A \bigg[ \dfrac{d}{dt}\bigg(\dfrac{mw}{|e|B}\bigg)\hat a +  \dfrac{mw}{|e|B}\dfrac{d\hat a}{dt}\bigg] \\
    &= e\dfrac{m}{|e|}\dfrac{mw}{|e|B}(\hat a\vec\nabla)\vec A \bigg[ \overbrace{\dfrac{\dot w B - {w\dot B}}{B^2}}^{\text{αργή μεταβολή του Β}}\hat a +  \dfrac{w}{B}\dot\xi\hat c\bigg] \\
    &=e\dfrac{m}{|e|}\dfrac{mw}{|e|B}\dfrac{\dot w}{B}(\hat a \vec\nabla)\vec A\cdot\hat a 
    + e\dfrac{m}{|e|}\dfrac{mw}{|e|B}\dfrac{w}{B}\dot\xi\cdot(\hat a \vec\nabla)\vec A\cdot\hat c\\
    &= + \dfrac{m^2}{e}\dfrac{w\dot w}{B^2}(\hat a \vec\nabla)\vec A\cdot\hat a 
    + \dfrac{m^2}{e}\dfrac{w^2}{B^2}\ \dot\xi\cdot(\hat a \vec\nabla)\vec A\cdot\hat c
\end{aligned}
$$

>Στην τελευταία έκφραση, οι συντελεστές έχουν μονάδες $[m^2Cs^{-1}]$ ενώ το γινόμενο διανυσμάτων έχει μονάδες $[kgC^{-1}s^{-1}]$ και άρα το γινόμενό τους έχει μονάδες $[kgm^2/s^2]$, δηλαδή ενέργειας, όπως θα έπρεπε.

Άρα η Λαγκραντζιανή γράφεται:

$$
\large\boxed{
    \begin{aligned}
    \mathcal{L} &= \bigg[e\vec A + m\upsilon_{||}\hat b\bigg]\dot{\vec X} 
    + \dfrac{m^2}{e}\dfrac{w^2}{B}\dot\xi 
    - e \bigg(\dfrac{mw}{|e|B}\bigg) \hat a \dfrac{\partial \vec A}{\partial t}\\
    &+ \dfrac{m^2}{e}\dfrac{w\dot w}{B^2}(\hat a \vec\nabla)\vec A\cdot\hat a 
    + \dfrac{m^2}{e}\dfrac{w^2}{B^2}\ \dot\xi\cdot(\hat a \vec\nabla)\vec A\cdot\hat c-\mathcal{H}
\end{aligned}
}
$$

---

Προσθέτουμε τώρα την (νέα) $\dfrac{dS}{dt}$, με $S = -\dfrac{m^2w^2}{2qB^2}(\hat a \vec\nabla)\vec A\cdot\hat a$

$$
\begin{aligned}
    \dfrac{dS}{dt} &= \dfrac{d}{dt}\bigg[ -\dfrac{m^2w^2}{2qB^2}(\hat a \vec\nabla)\vec A\cdot\hat a \bigg] =\\
    &=-\dfrac{d}{dt}\bigg(\dfrac{m^2w^2}{2qB^2}\bigg)(\hat a \vec\nabla)\vec A\cdot\hat a
    - \dfrac{m^2w^2}{2qB^2}\dfrac{d\hat a}{dt}\vec\nabla\vec A\cdot\hat a
    - \overbrace{\dfrac{m^2w^2}{2qB^2}\hat a\dfrac{d(\vec\nabla\vec A)}{dt}\hat a}^{\approx\mathcal{O}(\epsilon^2)}
    - \dfrac{m^2w^2}{2qB^2}(\hat a \vec\nabla)\vec A \dfrac{d\hat a}{dt} = \\
    &=-\bigg(\dfrac{m^2}{2q}\bigg)\dfrac{2w\dot w B^2 - \overbrace{2w^2B\dot B}^{\approx\mathcal{O}(\epsilon^2)}}{B^4}(\hat a \vec\nabla)\vec A\cdot\hat a
    -\dfrac{m^2w^2}{2qB^2}\ \dot\xi\ \hat c\ \vec\nabla\vec A\cdot\hat a
    -\overbrace{\dfrac{m^2w^2}{2qB^2}(\hat a \vec\nabla)\vec A (\dot\xi\ \hat c)}^{\text{προσθαφαιρούμε }\times2}=\\ \\
    &=-\dfrac{m^2}{e}\dfrac{w\dot w}{B^2}(\hat a \vec\nabla)\vec A\cdot\hat a
    -\dfrac{m^2}{e}\dfrac{w^2}{2B^2}\dot\xi(\hat c\vec\nabla)\vec A\cdot\hat a
    +\dfrac{m^2}{e}\dfrac{w^2}{2B^2}\ \dot\xi (\hat a\vec\nabla) \vec A \hat c
    -2 \dfrac{m^2}{e}\dfrac{w^2}{2B^2}\ \dot\xi (\hat a\vec\nabla) \vec A \hat c=\\
\end{aligned} 
$$

Εξ'ορισμού, έχουμε $\partial_xA_y - \partial_yA_x = B_z + $ μεταθέσεις. Αν αντί για $(\hat x, \hat y,\hat z)$ χρησιμοποιήσουμε $(\hat c, \hat a, \hat b)$ και προσθέσουμε όλες τις μεταθέσεις διανυσματικά, προκύπτει η σχέση $(\hat x \vec\nabla)\vec A\cdot\hat a - (\hat a \vec\nabla)\vec A\cdot\hat c  = B$, επομένως:

$$
\begin{aligned}\dfrac{dS}{dt} 
    &=-\dfrac{m^2}{e}\dfrac{w\dot w}{B^2}(\hat a \vec\nabla)\vec A\cdot\hat a
    -\dfrac{m^2}{e}\dfrac{w^2}{2B^2}\dot\xi\bigg[\overbrace{ (\hat c\vec\nabla)\vec A\cdot\hat a - (\hat a\vec\nabla) \vec A \hat c }^B\bigg]
    -2 \dfrac{m^2}{e}\dfrac{w^2}{2B^2}\ \dot\xi (\hat a\vec\nabla) \vec A \hat c
    \Rightarrow\\ \\
    \dfrac{dS}{dt} 
    &=\overbrace{-\dfrac{m^2}{e}\dfrac{w\dot w}{B^2}(\hat a \vec\nabla)\vec A\cdot\hat a}^{\text{θα φύγει}}
    -\dfrac{m^2}{e}\dfrac{w^2}{2B}\dot\xi
    -\overbrace{\dfrac{m^2}{e}\dfrac{w^2}{B^2}\ \dot\xi (\hat a\vec\nabla) \vec A \hat c}^{\text{θα φύγει}}
\end{aligned} 
$$

> Και εδώ, και οι τρεις όροι της $dS/dt$ έχουν μονάδες ενέργειας όπως θα έπρεπε.

---

Αμέσως υπολογίζεται η νέα Λαγκραντζιανή: $\mathcal{L} = \mathcal{L}_{\text{παλιά}} + \dfrac{dS}{dt}$

$$
\begin{aligned}
    \mathcal{L} &= \bigg[e\vec A + m\upsilon_{||}\hat b\bigg]\dot{\vec X} 
    + \dfrac{m^2}{e}\dfrac{w^2}{B}\dot\xi 
    - e \bigg(\dfrac{mw}{qB}\bigg) \hat a \dfrac{\partial \vec A}{\partial t}
    -\dfrac{m^2}{e}\dfrac{w^2}{2B}\dot\xi
    -\mathcal{H} \\
    &=\bigg[e\vec A + m\upsilon_{||}\hat b\bigg]\dot{\vec X} 
    - \dfrac{m^2}{e}\dfrac{w^2}{2B}\dot\xi 
    - e \bigg(\dfrac{mw}{qB}\bigg) \hat a \dfrac{\partial \vec A}{\partial t}
    -\mathcal{H}
\end{aligned}
$$

Τέλος, παίρνοντας τον μέσο όρο της κυκλοτρονικής κίνησης, απαλοίφεται ο όρος $\partial \vec A/dt$, οπότε προκύπτει η δεύτερης τάξης Λαγκραντζιανή του γυροκέντρου:

$$
\large\boxed{
    \mathcal{L}(\vec X,\upsilon_{||},w,\dot\xi,t)
    =\bigg[e\vec A + m\upsilon_{||}\hat b\bigg]\dot{\vec X} 
    - \dfrac{m^2}{e}\dfrac{w^2}{2B}\dot\xi 
    -\mathcal{H}(\vec X,\upsilon_{||},w,t)
}
$$

---

Για να βρούμε τη Χαμιλτονιανή, ξεκινάμε από τον αρχικό τύπο:

$$
\begin{aligned}
    \mathcal{H} 
    &= \dfrac{1}{2}m\vec u^2 + e\Phi(\vec x,t)
    =\dfrac{1}{2}m(\upsilon_{||}^2 + w^2) + e\Phi(\vec x,t)\\
    &=\dfrac{1}{2}m\upsilon_{||}^2 + \dfrac{1}{2}mw^2 + e\Phi(\vec x,t)
    =\dfrac{1}{2}m\upsilon_{||}^2 + \overbrace{m\dfrac{w^2}{2B}}^\mu B + e\Phi(\vec x,t)
\end{aligned}
$$

Και τελικά η Χαμιλτονιανή γράφεται ως:

$$
\large\boxed{
    \mathcal{H}(\vec X,\upsilon_{||},w,t)
    =\dfrac{1}{2}m\upsilon_{||}^2 + \mu B + e\Phi(\vec X,t)
    \quad,\quad 
    \mu = m\dfrac{w^2}{2B} = \text{σταθ}
}
$$

## Λαγκραντζιανός και Χαμιλτονιανός Φορμαλισμός

Ξαναγράφουμε την Λαγκραντζιανή και Χαμιλτονιανή του σωματιδίου χρησιμοποιώντας:

1. Το $\rho_{||}$ είναι η "παράλληλη γυροακτίνα", και είναι κανονικοποιημένη ως προς την μεγάλη ακτίνα $R$. Έχει μονάδες $[m^2T^{-1}s^{-1}]$ και ορίζεται ώς  
 $\quad \rho_{||} = \dfrac{R\upsilon_{||}}{B},\quad$ και άρα $\quad m\dfrac{\rho_{||}}{R}\vec B = m\upsilon_{||}\hat b, \quad$ και $\quad \upsilon_{||}^2 = \bigg(\dfrac{\rho_{||}}{R}\bigg)^2B^2$

1. $\quad \mu = m\dfrac{w^2}{2B},\quad$ και άρα $\quad \dfrac{m^2}{e}\dfrac{w^2}{2B}\dot\xi = \dfrac{m}{e}\mu\dot\xi$

2. $\quad \dot{\vec X} = \vec \upsilon$

Επομένως:

$$
\large\boxed{
    \mathcal{L}
    =\bigg[e\vec A + m\dfrac{\rho_{||}}{R}\vec B\bigg]\vec \upsilon 
    - \dfrac{m}{e}\mu\dot\xi
    -\mathcal{H}
}
$$

<center>
και
</center>

$$
\large\boxed{
    \mathcal{H}
    =\dfrac{1}{2}m\bigg(\dfrac{\rho_{||}}{R}\bigg)^2B^2 + \mu B + e\Phi
}
$$

Ξαναγράφουμε την ανταλλοίωτη μορφή του μαγνητικού πεδίου:

$$
\qquad \vec B = \nabla\psi \times \nabla\theta - \nabla\psi_p \times \nabla\zeta
= \nabla\times\overbrace{(\psi\nabla\theta - \psi_p\nabla\zeta)}^{\vec A} = \nabla\vec A
$$

και αντικαθιστούμε στην Λαγκραντζιανή:

$$
\begin{aligned} \mathcal{L}
    &=\bigg[e\overbrace{(\psi\nabla\theta - \psi_p\nabla\zeta)}^{\vec A} + m\mu_0\dfrac{\rho_{||}}{R}\overbrace{(g(\psi)\nabla\zeta + I(\psi)\nabla\theta + \delta(\psi,\theta)\nabla\psi)}^{\vec B}\bigg]\vec\upsilon  - \dfrac{m}{e}\mu\dot\xi-\mathcal{H}
\end{aligned}
$$

Ακόμα, θα ισχύει

$$
\vec\upsilon = \upsilon_\psi\nabla\psi + \upsilon_\theta\nabla\theta + \upsilon_\zeta\nabla\zeta
=\dot\psi\nabla\psi + \dot\theta\nabla\theta + \dot\zeta\nabla\zeta
$$

Υπολογίζουμε τους δύο όρους της επιμεριστικής ξεχωριστά:

$$
\begin{aligned}
    \bullet\quad e\vec A \cdot \vec \upsilon
    &= e\bigg[\psi\nabla\theta - \psi_p\nabla\zeta\bigg]\cdot\bigg[\dot\psi\nabla\psi + \dot\theta\nabla\theta + \dot\zeta\nabla\zeta\bigg]
    =e\psi\dot\theta - e\psi_p\dot\zeta
\end{aligned}
$$

$$
\begin{aligned}
    \bullet\quad m\mu_0\dfrac{\rho_{||}}{R}\vec B \cdot \vec\upsilon
    &= m\mu_0\dfrac{\rho_{||}}{R}\bigg[g\nabla\zeta + I\nabla\theta + \delta\nabla\psi\bigg]\cdot
    \bigg[\dot\psi\nabla\psi + \dot\theta\nabla\theta + \dot\zeta\nabla\zeta\bigg] \\
    &=m\mu_0\dfrac{\rho_{||}}{R}\bigg[ g\dot\zeta + I\dot\theta + \delta\dot\psi\bigg]
\end{aligned}
$$

Από τον ορισμό του $\psi_p$ έχουμε:

$$
\begin{aligned} \dot\psi 
= \dfrac{d\psi}{dt} = \dfrac{d\psi}{d\psi_p}\dfrac{d\psi_p}{dt}
= q\dot\psi_p
\end{aligned}
$$

Άρα:

$$
\begin{aligned}
    \bullet\quad m\mu_0\dfrac{\rho_{||}}{R}\vec B \cdot \vec\upsilon
    =m\mu_0\dfrac{\rho_{||}}{R}g\dot\zeta + m\mu_0\dfrac{\rho_{||}}{R}I\dot\theta + m\mu_0\dfrac{\rho_{||}}{R}\delta q\dot\psi_p
\end{aligned}
$$
 
Και επομένως η Λαγκραντζιανή θα γράφεται:

$$
\mathcal{L}= e\psi\dot\theta - e\psi_p\dot\zeta + m\mu_0\dfrac{\rho_{||}}{R}g\dot\zeta + m\mu_0\dfrac{\rho_{||}}{R}I\dot\theta + m\mu_0\delta\dfrac{\rho_{||}}{R} - \dfrac{m}{e}\mu\dot\xi-\mathcal{H}
\Rightarrow 
$$

$$
\large\boxed{
    \mathcal{L}
    =(e\psi + m\mu_0\dfrac{\rho_{||}}{R}I)\dot\theta + (m\mu_0\dfrac{\rho_{||}}{R}g - e\psi_p)\dot\zeta
    - \dfrac{m}{e}\mu\dot\xi + m\mu_0\dfrac{\rho_{||}}{R}\delta q \dot\psi_p - \mathcal{H}
}
$$

> Οι μαγνητικές ροές $\psi$ και $\psi_p$ έιναι σε $[Tm^2]$, το $\mu_0$ σε $[ΝΑ^{-2}]$, τα ρεύματα $g,I$ σε $[Α]$, και άρα όλοι οι όροι έχουν μονάδες ενέργειας, όπως θα έπρεπε.

Τώρα έχουμε φέρει τη Λαγκραντζιανή στη μορφή

$$
\mathcal{L} (\vec q, \dot{\vec q}) = \mathcal{L} (\rho_{||}, \psi_p, \theta, \zeta,\quad \dot\rho_{||}, \dot\psi_p, \dot\theta, \dot\zeta) 
$$ 

και μπορούμε να βρούμε τις εξισώσεις κίνησης μέσω των εξισώσεων Euler-Lagrange.

---

### Εξισώσεις Κίνησης

Αγνοούμε τον όρο $\delta$, καθώς επηρεάζει μόνο την χρονική κλίμακα και όχι την κίνηση. Ακόμα, θεωρούμε ότι τα ρεύματα $g,I$ είναι σταθερά και δεν εξαρτώνται από καμία συντεταγμένη. Άρα η Λαγκραντζιανή που θα χρησιμοποιούμε απο εδώ και πέρα είναι η:

$$
\large\boxed{
    \mathcal{L}
    =(e\psi + \dfrac{m\mu_0\rho_{||}I}{R})\dot\theta + (\dfrac{m\mu_0\rho_{||}g}{R} - e\psi_p)\dot\zeta
    - \dfrac{m}{e}\mu\dot\xi  - \mathcal{H}
}
$$

Εξισώσεις Euler-Lagrange:

$$
\begin{aligned}
    \dfrac{d}{dt}\dfrac{\partial\mathcal{L}}{\partial \dot\rho_{||}} = \dfrac{\partial\mathcal{L}}{\partial \rho_{||}} \quad&\Rightarrow\quad
    \dfrac{m\mu_0g}{R}\dot\zeta + \dfrac{m\mu_0I}{R}\dot\theta  =\partial_{\rho_{||}}\mathcal{H} 
    \\[15pt]
    \dfrac{d}{dt}\dfrac{\partial\mathcal{L}}{\partial \dot\psi_p} = \dfrac{\partial\mathcal{L}}{\partial \psi_p} \quad&\stackrel{\dot\psi = q\dot\psi_p}{\Rightarrow}\quad
    eq\dot\theta - e\dot\zeta = \partial_{\psi_p}\mathcal H
    \\[15pt]
    \dfrac{d}{dt}\dfrac{\partial\mathcal{L}}{\partial \dot\theta} = \dfrac{\partial\mathcal{L}}{\partial \theta} \quad&\stackrel{\dot\psi = q\dot\psi_p}{\Rightarrow}\quad 
    -eq\dot\psi_p - \dfrac{m\mu_0\dot\rho_{||}I}{R} = \partial_\theta\mathcal{H}
    \\[15pt]
    \dfrac{d}{dt}\dfrac{\partial\mathcal{L}}{\partial \dot\zeta} = \dfrac{\partial\mathcal{L}}{\partial \zeta} \quad&\Rightarrow\quad
    \dfrac{m\mu_0\dot\rho_{||}g}{R}- + e\dot\psi_p = \partial_\zeta\mathcal{H}
\end{aligned}
$$

Η σε μορφή πινάκων:

$$
\begin{bmatrix}
0 & 0 & \dfrac{m\mu_0 I}{R} & \dfrac{m\mu_0 g}{R} \\[10pt]
0 & 0 & eq & -e \\[10pt]
-\dfrac{m\mu_0 I}{R} & -eq & 0 & 0 \\[10pt]
-\dfrac{m\mu_0 g}{R} & e & 0 & 0
\end{bmatrix}
\cdot
\begin{bmatrix}
\dot\rho_{||} \\
\dot\psi_p \\
\dot\theta \\
\dot\zeta
\end{bmatrix}
=
\begin{bmatrix}
    \partial_{\rho_{||}}\mathcal{H} \\
    \partial_{\psi_p}\mathcal{H} \\
    \partial_\theta\mathcal{H} \\
    \partial_\zeta\mathcal{H}
\end{bmatrix}
$$

Έχουμε έναν πίνακα της μορφής $Μ = \begin{bmatrix} 0 & B \\ C & 0 \end{bmatrix}$. O αντίστροφος ενός τέτοιου πίνακα είναι ο $Μ^{-1} = \begin{bmatrix} 0 & C^{-1} \\ B^{-1} & 0 \end{bmatrix}$ και η διακρίνουσά του $|Μ|$ θα είναι:

$$
\begin{aligned}
|Μ| &= |A||B| = (-\dfrac{me\mu_0 I}{R} - \dfrac{me\mu_0 gq}{R})(-\dfrac{me\mu_0I}{R} - \dfrac{me\mu_0 gq}{R}) \\[10pt]
&= \big[(\dfrac{me\mu_0}{R})(I+gq)\big]^2   = D^2 \neq 0
\end{aligned}
$$

Άρα μπορούμε να αντιστρέψουμε το παραπάνω σύστημα των εξισώσεων κίνησης:

$$
\begin{bmatrix}
    \dot\rho_{||} \\
    \dot\psi_p \\
    \dot\theta \\
    \dot\zeta
\end{bmatrix}
= \dfrac{1}{D}
\begin{bmatrix}
    0 & 0 & -e & -eq \\[10pt]
    0 & 0 & \dfrac{m\mu_0 g}{R} & \dfrac{m\mu_0 I}{R}\\[10pt]
    e & \dfrac{m\mu_0 g}{R} & 0 & 0 \\[10pt]
    eq & \dfrac{m\mu_0 I}{R} & 0 & 0
\end{bmatrix}
\cdot
\begin{bmatrix}
    \partial_{\rho_{||}}\mathcal{H} \\
    \partial_{\psi_p}\mathcal{H} \\
    \partial_\theta\mathcal{H} \\
    \partial_\zeta\mathcal{H}
\end{bmatrix}
$$

> Οι ορίζουσες που προκύπτουν από την αντιστροφή των πινάκων $A,B$ είναι ίσοι με $-D$, για αυτό εμφανίζεται μόνο ο όρος $1/D$ και τα πρόσημα είναι αντεστραμμένα.

Η Χαμιλτονιανή βρήκαμε πως είναι: 

$$
\large\boxed{
    \mathcal{H}
    =\dfrac{1}{2}m\bigg(\dfrac{\rho_{||}}{R}\bigg)^2B^2 + \mu B + e\Phi
}
$$
Υπολογίζουμε λοιπόν τις ποσότητες:

$$
\begin{aligned}
    \partial_{\rho_{||}}\mathcal{H} &= m\bigg(\dfrac{\rho_{||}}{R}\bigg)^2B^2 \\[10pt]
    \partial_{\psi_p}\mathcal{H} 
    &= \bigg(\mu + m\bigg(\dfrac{\rho_{||}}{R}\bigg)^2B\bigg)\dfrac{\partial B}{\partial \psi_p} 
    + e\dfrac{\partial \Phi}{\partial \psi_p} \\[12pt]
    \partial_\theta\mathcal{H}
    &= \bigg(\mu + m\bigg(\dfrac{\rho_{||}}{R}\bigg)^2B\bigg)\dfrac{\partial B}{\partial \theta}
    + e\dfrac{\partial \Phi}{\partial \theta} \\[12pt]
    \partial_\zeta\mathcal{H} &= 0
\end{aligned}
$$

Με απλή αντικατάσταση προκύπτουν τώρα οι εξισώσεις κίνησης:

$$
\large\bullet\quad \dot\theta
    = em\dfrac{\rho_{||}B^2}{D}
    + m\mu_0\dfrac{g}{RD}\bigg[\bigg(\mu + m\bigg(\dfrac{\rho_{||}}{R}\bigg)^2B\bigg)\dfrac{\partial B}{\partial \psi_p} 
    + e\dfrac{\partial \Phi}{\partial \psi_p}\bigg]
    \\[20pt]
\large\bullet\quad \dot \psi_p
    = -m\mu_0\dfrac{g}{RD}\bigg[\bigg(\mu + m\bigg(\dfrac{\rho_{||}}{R}\bigg)^2B\bigg)\dfrac{\partial B}{\partial \theta}
    + e\dfrac{\partial \Phi}{\partial \theta}\bigg]
    \\[20pt]
\large\bullet\quad \dot\rho_{||}
    = -e\dfrac{1}{RD}\bigg[\bigg(\mu + m\bigg(\dfrac{\rho_{||}}{R}\bigg)^2B\bigg)\dfrac{\partial B}{\partial \theta} + e\dfrac{\partial \Phi}{\partial \theta}\bigg]
    \\[20pt]
\large\bullet\quad \dot\zeta
    = eqm\dfrac{\rho_{||}B^2}{D} 
    - m\mu_0\dfrac{I}{RD}\bigg[\bigg(\mu + m\bigg(\dfrac{\rho_{||}}{R}\bigg)^2B\bigg)\dfrac{\partial B}{\partial \psi_p} 
    + e\dfrac{\partial \Phi}{\partial \psi_p}\bigg]
    \\[20pt]
\text{όπου}\quad D = \dfrac{me\mu_0}{R}(I+gq)
$$
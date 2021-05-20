.. default-role:: math


Unit conversion tool
--------------------

This is an automated unit conversion tool. 
Fill in the 4 fundamental quantities, `\hat{c}`, (numerical speed of light), `\hat{\mathcal{R}}` (skin depth resolution), `\hat{n}_{ppc}` (particles per cell per species), `\sigma` (magnetization), and the resulting numerical and physical values are automatically updated.

Add also a length scale, `\Delta x/1 \mathrm{cm}`, to get the resulting physical values.

----


.. raw:: html

  <div>

    <div class="row">
      <div class="column">
        <div id="c-input" class="independent-value setvalue-group col-xs-12">
          <div class="value-input">
            <p> <span class="math notranslate nohighlight">\( \hat{c}= \)</span> <input class="value input-small" value="0.45"></p>
          </div>
        </div>
      </div>

      <div class="column">
        <div id="comp-input" class="independent-value setvalue-group col-xs-12">
          <div class="value-input">
            <p> <span class="math notranslate nohighlight">\( \hat{\mathcal{R}}= \)</span> <input class="value input-small" value="5"></p>

          </div>
        </div>
      </div>

      <div class="column">
        <div id="nppc-input" class="independent-value setvalue-group col-xs-12">
          <div class="value-input">
            <p> <span class="math notranslate nohighlight">\( \hat{n}_{\mathrm{ppc}} = \)</span> <input class="value input-small" value="1"></p>

          </div>
        </div>
      </div>

      <div class="column">
        <div id="sigma-input" class="independent-value setvalue-group col-xs-12">
          <div class="value-input">
            <p> <span class="math notranslate nohighlight">\( \sigma= \)</span> <input class="value input-small" value="10"></p>
          </div>
        </div>
      </div>
    </div>

    <div class="row">
      <div class="column" style="width:100%;float:left;">
        <div id="dx-input" class="independent-value setvalue-group col-xs-12">
          <div class="value-input">
            <p> <span class="math notranslate nohighlight">\( \Delta x/(1~\mathrm{cm}) = \)</span> <input class="value input-small" value="1"></p>
          </div>
        </div>
      </div>
    </div>

  </div>


-----

Numerical simulation variables
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


.. raw:: html

    <!-- OUTPUT -->
    <div class="row">
      <div class="column" style="width:33%">
        <div id="omegap-output" class="dependent-value getvalue-group col-xs-12">
          <div class="value-output">
            
            <p>
             <span class="math notranslate nohighlight">\( \omega_{pe}= \)</span>
             <span class="value" style="color:blue;"></span> 
             <span class="math notranslate nohighlight">\( \Delta t^{-1} \)</span>
             </p>

          </div>
        </div>
      </div>

      <div class="column" style="width:33%">
        <div id="omegaB-output" class="dependent-value getvalue-group col-xs-12">
          <div class="value-output">
        
             <p>
             <span class="math notranslate nohighlight">\( \omega_B= \)</span>
             <span class="value" style="color:blue;"></span> 
             <span class="math notranslate nohighlight">\( (\gamma\beta)^{-1} \Delta t^{-1}\)</span>
             </p>
             
          </div>
        </div>
      </div>

      <div class="column" style="width:33%">
        <div id="gyro-output" class="dependent-value getvalue-group col-xs-12">
          <div class="value-output">
            
             <p>
             <span class="math notranslate nohighlight">\( r_L= \)</span>
             <span class="value" style="color:blue;"></span> 
             <span class="math notranslate nohighlight">\( \gamma\beta \Delta x\)</span>
             </p>
             
          </div>
        </div>
      </div>
    </div>

    <div style="width:100%;float:left;"> </div>

    <div class="row">

      <div class="column" style="width:33%;float:left;">
        <div id="valf-output" class="dependent-value getvalue-group col-xs-12">
          <div class="value-output">
            
             <p>
             <span class="math notranslate nohighlight">\( v_A= \)</span>
             <span class="value" style="color:blue;"></span> 
             <span class="math notranslate nohighlight">\( c \)</span>
             </p>
             
          </div>
        </div>
      </div>

      <div class="column" style="width:33%;float:left;">
        <div id="binit-output" class="dependent-value getvalue-group col-xs-12">
          <div class="value-output">
            
             <p>
             <span class="math notranslate nohighlight">\( B_0 = \)</span>
             <span class="value" style="color:blue;"></span> 
             </p>
             
          </div>
        </div>
      </div>

      <div class="column" style="width:33%:float:left;">
        <div id="qe-output" class="dependent-value getvalue-group col-xs-12">
          <div class="value-output">
            
             <p>
             <span class="math notranslate nohighlight">\( | q_e | = m_e = \)</span>
             <span class="value" style="color:blue;"></span> 
             </p>
             
          </div>
        </div>
      </div>
    </div>
  </div>
  

  <div style="width:100%;float:left;"> </div>



------

.

Physical scales 
^^^^^^^^^^^^^^^

.. raw:: html

    <!-- OUTPUT -->

    <div style="width:100%;float:left;"> </div>
    <div class="row">
      <div class="column" style="width:50%">
        <div id="bphys-output" class="dependent-value getvalue-group col-xs-12">
          <div class="value-output">
            
            <p>
             <span class="math notranslate nohighlight">\( B = \)</span>
             <span class="value" style="color:red;"></span> 
             <span class="math notranslate nohighlight">\( \hat{B} ~\mathrm{Gauss} \)</span>
             </p>

          </div>
        </div>
      </div>

      <div class="column" style="width:50%">
        <div id="ephys-output" class="dependent-value getvalue-group col-xs-12">
          <div class="value-output">
            
            <p>
             <span class="math notranslate nohighlight">\( E = \)</span>
             <span class="value" style="color:red;"></span> 
             <span class="math notranslate nohighlight">\( \hat{E} ~\mathrm{statvolt}~\mathrm{cm}^{-1} \)</span>
             </p>

          </div>
        </div>
      </div>
    </div>

  <div style="width:100%;float:left;"> </div>

    <div class="row">
      <div class="column" style="width:50%">
        <div id="jphys-output" class="dependent-value getvalue-group col-xs-12">
          <div class="value-output">
            
            <p>
             <span class="math notranslate nohighlight">\( J = \)</span>
             <span class="value" style="color:red;"></span> 
             <span class="math notranslate nohighlight">\( \hat{J} ~\mathrm{statcoul}~\mathrm{s}^{-1} \)</span>
             </p>

          </div>
        </div>
      </div>

      <div class="column" style="width:50%">
        <div id="qphys-output" class="dependent-value getvalue-group col-xs-12">
          <div class="value-output">
            
            <p>
             <span class="math notranslate nohighlight">\( q_e = \)</span>
             <span class="value" style="color:red;"></span> 
             <span class="math notranslate nohighlight">\( ~\mathrm{statcoul} \)</span>
             </p>

          </div>
        </div>
      </div>
    </div>

  <div style="width:100%;float:left;"> </div>

    <div class="row">

      <div class="column" style="width:33%">
        <div id="skphys-output" class="dependent-value getvalue-group col-xs-12">
          <div class="value-output">
            
            <p>
             <span class="math notranslate nohighlight">\( \frac{c}{\omega_{pe}} = \)</span>
             <span class="value" style="color:red;"></span> 
             <span class="math notranslate nohighlight">\( ~\mathrm{cm} \)</span>
             </p>

          </div>
        </div>
      </div>

      <div class="column" style="width:33%">
        <div id="omphys-output" class="dependent-value getvalue-group col-xs-12">
          <div class="value-output">
            
            <p>
             <span class="math notranslate nohighlight">\( \omega_{pe} = \)</span>
             <span class="value" style="color:red;"></span> 
             <span class="math notranslate nohighlight">\( ~\mathrm{s}^{-1} \)</span>
             </p>

          </div>
        </div>
      </div>

      <div class="column" style="width:33%">
        <div id="nphys-output" class="dependent-value getvalue-group col-xs-12">
          <div class="value-output">
            
            <p>
             <span class="math notranslate nohighlight">\( n_e = \)</span>
             <span class="value" style="color:red;"></span> 
             <span class="math notranslate nohighlight">\( ~\mathrm{cm}^{-3} \)</span>
             </p>

          </div>
        </div>
      </div>
    </div>

  <div style="width:100%;float:left;"> </div>

.. raw:: html

  <script>
    window.onload = function() {

      updateDependents();


      let independent_values = document.getElementsByClassName("independent-value");
      [].forEach.call(independent_values, function(value) {
        let value_inputs = value.getElementsByClassName("value");
        [].forEach.call(value_inputs, function(value_input) {
          value_input.addEventListener("input", updateDependents, false);
        });
      });

      function precise(x) {
        if (isFinite(x)) {
          return Number.parseFloat(x).toPrecision(4);
        } else {
          return '&#8734';
        }
      }

      function expo_precise(x) {
        if (isFinite(x)) {
          return Number.parseFloat(x).toExponential(4);
        } else {
          return '&#8734';
        }
      }


      function updateDependents() {
        let cfl = 1.0*document.getElementById("c-input").getElementsByClassName("value")[0].value;
        let comp = 1.0*document.getElementById("comp-input").getElementsByClassName("value")[0].value;
        let sigma = 1.0*document.getElementById("sigma-input").getElementsByClassName("value")[0].value;
        let ppc = 1.0*document.getElementById("nppc-input").getElementsByClassName("value")[0].value;
        let dx = 1.0*document.getElementById("dx-input").getElementsByClassName("value")[0].value;

        let dt = dx/3e10;

        let qe = cfl ** 2 / (2*ppc * comp ** 2);
        let me = qe;

        let b_norm = Math.sqrt( 2*ppc*cfl*cfl*me*sigma );

        let omega_P = cfl/comp;
        let omega_B = Math.sqrt(sigma)*omega_P;

        let v_A = Math.sqrt(sigma/(sigma+1));


        {
          // update omega_p
          let omegap_el = document.getElementById("omegap-output");
          omegap_el.getElementsByClassName("value")[0].innerHTML = precise(omega_P);
        }

        {
          // update omegaB
          let omegap_B = document.getElementById("omegaB-output");
          omegap_B.getElementsByClassName("value")[0].innerHTML = precise(omega_B);
        }

        {
          // update gyro
          let gyro_el = document.getElementById("gyro-output");
          gyro_el.getElementsByClassName("value")[0].innerHTML = precise(comp / Math.sqrt(sigma));
        }

        {
          // update alfven
          let valf = document.getElementById("valf-output");
          valf.getElementsByClassName("value")[0].innerHTML = precise( v_A  );
        }

        {
          // update B0
          let binit = document.getElementById("binit-output");
          binit.getElementsByClassName("value")[0].innerHTML = precise( b_norm  );
        }

        {
          // update qe
          let qev = document.getElementById("qe-output");
          qev.getElementsByClassName("value")[0].innerHTML = precise( qe  );
        }


        {
          // update b
          let bphys = document.getElementById("bphys-output");
          bphys.getElementsByClassName("value")[0].innerHTML = expo_precise( 1.705e3*(cfl**2)/dx  );
        }

        {
          // update e
          let ephys = document.getElementById("ephys-output");
          ephys.getElementsByClassName("value")[0].innerHTML = expo_precise( 1.705e3*(cfl**2)/dx  );
        }

        {
          // update j
          let jphys = document.getElementById("jphys-output");
          jphys.getElementsByClassName("value")[0].innerHTML = expo_precise( 4.056e12*(cfl**3)/dx**2  );
        }

        {
          // update q
          let qphys = document.getElementById("qphys-output");
          qphys.getElementsByClassName("value")[0].innerHTML = expo_precise( 1.356e2*qe*(cfl**2)*dx  );
        }

        {
          // update sk
          let skphys = document.getElementById("skphys-output");
          skphys.getElementsByClassName("value")[0].innerHTML = expo_precise( comp*dx  );
        }

        {
          // update omegap
          let omphys = document.getElementById("omphys-output");
          omphys.getElementsByClassName("value")[0].innerHTML = expo_precise( omega_P*dt  );
        }

        {
          // update numdens
          let nphys = document.getElementById("nphys-output");

          let bnorm = 1.705e3*(cfl**2)/dx;

          //b2/4*pi*n*me*c^2 = sigma
          //n = 4.0*3.14*((3.0e10)**2)*9.1e-28/b^2
          //n = 1.028e-5/b^2
          let nphysv = (bnorm**2)/(1.028e-5*sigma);

          nphys.getElementsByClassName("value")[0].innerHTML = expo_precise( nphysv );
        }

      }

    };

  </script>







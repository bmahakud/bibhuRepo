import React,{useState} from 'react';



import classes from './App.module.css';

function App() {


  const [homepage, setHomePage]= useState(true);
  


  function showHomePage(){

    setHomePage(true);

  }


  function showAboutPage(){

    setHomePage(false);

  }




  return (
    <div className={classes.app}>

      <div className={classes.header} >
        <div style={{borderStyle:'solid',color:"blue",marginLeft:"20px"}} 
             onClick={showHomePage}
             >
          home
        </div>

        <div style={{borderStyle:'solid',color:'blue',marginLeft:"20px" }}
           onClick={showAboutPage}
        >
          about

        </div>           

        
      </div>
    

      <div className={classes.body}>


        { homepage &&
        <div className={classes.home}>
    
           This is home page  
         
        </div>

        }


        { !homepage &&
        <div className={classes.about}>

          This is about page
        </div>

        }


      </div>





    </div>
  );
}

export default App;
